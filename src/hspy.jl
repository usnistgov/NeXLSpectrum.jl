using HDF5: attributes, h5open, ishdf5, File
using TimeZones: TimeZone, ZonedDateTime, localzone

"""
    readhspy(filename::String, name=nothing)

Read a HyperSpectrum from a HyperSpy-style HDF5 file.  The `name` argument allows the user to
optionally specify an experiment to load.
"""
function readhspy(filename::String, name=nothing)
    function read_attributes(fid)
        cs = attributes(fid)
        Dict(map(k->k=>read(cs,k),keys(cs))...)
    end
    function assign(props, prop, val) 
        if !ismissing(val) props[prop]=val end
    end
    function read_stage(stage_attr)
        Dict(
            :X => get(stage_attr, "x", 0.0)*0.1,
            :Y => get(stage_attr, "y", 0.0)*0.1,
            :Z => get(stage_attr, "z", 0.0)*0.1,
            :R => deg2rad(get(stage_attr, "rotation", 0.0)),
            :T => deg2rad(get(stage_attr, "tilt_alpha", 0.0)),
            :B => deg2rad(get(stage_attr, "tilt_beta", 0.0)),
        )
    end
    function read_eds(props, det_attr)
        assign(props, :Azimuthal, get(det_attr, "azimuth_angle", missing)*deg2rad(1.0))
        assign(props, :Elevation, get(det_attr, "elevation_angle", missing)*deg2rad(1.0))
        assign(props, :HSResolution, get(det_attr, "energy_resolution_MnKa", missing))
        assign(props, :LiveTime, get(det_attr, "live_time", missing))
        assign(props, :RealTime, get(det_attr, "real_time", missing))
    end
    function read_datetime(general)
        tz = get(general, "time_zone", missing)
        res = nothing
        try
            dt = get(general, "date", "1.1.1900")
            res = DateTime(dt, dateformat"dd.mm.yyyy")
        catch
            @warn "Error parsing date: $dt not of format dd.mm.yyyy"
            res = Date("1.1.1970", dateformat("dd.mm.yy"))
        end
        try
            tm = get(general, "time", missing)
            if !ismissing(tm)
                res += Time(tm, dateformat"HH:MM:SS")
            end
        catch
            @warn "Error parsing time: $tm not of format HH:MM:SS"
        end
        try
            if !ismissing(tz)
                res = ZonedDateTime(res, TimeZone(tz))
            end
        catch
            @warn "Error parsing timezone: $tz is not a standard timezone."
        end
        res
    end
    #= 
        ├── Acquisition_instrument
        │   ├── SEM
        │   │   ├── Detector
        │   │   │   ├── detector_type
        │   │   │   └── EDS
        │   │   │       ├── azimuth_angle (º)
        │   │   │       ├── elevation_angle (º)
        │   │   │       ├── energy_resolution_MnKa (eV)
        │   │   │       ├── live_time (s)
        │   │   │       └── real_time (s)
        │   │   ├── beam_current (nA)
        │   │   ├── beam_energy (keV)
        │   │   ├── probe_area (nm²)
        │   │   ├── convergence_angle (mrad)
        │   │   ├── magnification
        │   │   ├── microscope
        │   │   ├── Stage
        │   │   │   ├── rotation (º)
        │   │   │   ├── tilt_alpha (º)
        │   │   │   ├── tilt_beta (º)
        │   │   │   ├── x (mm)
        │   │   │   ├── y (mm)
        │   │   │   └── z (mm)
        │   │   └── working_distance (mm)
        │   └── TEM
        │       ├── Detector
        │       │   ├── EDS
        │       │   │   ├── azimuth_angle (º)
        │       │   │   ├── elevation_angle (º)
        │       │   │   ├── energy_resolution_MnKa (eV)
        │       │   │   ├── live_time (s)
        │       │   │   └── real_time (s)
        │       │   └── EELS
        │       │       ├── aperture (mm)
        │       │       ├── collection_angle (mrad)
        │       │       ├── dwell_time (s)
        │       │       ├── exposure (s)
        │       │       ├── frame_number
        │       │       └── spectrometer
        │       ├── Biprism
        │       │   ├── azimuth_angle (º)
        │       │   ├── position
        │       │   └── voltage (V)
        │       ├── acquisition_mode
        │       ├── beam_current (nA)
        │       ├── beam_energy (keV)
        │       ├── probe_area (nm²)
        │       ├── camera_length (mm)
        │       ├── convergence_angle (mrad)
        │       ├── magnification
        │       ├── microscope
        │       └── Stage
        │           ├── rotation (º)
        │           ├── tilt_alpha (º)
        │           ├── tilt_beta (º)
        │           ├── x (mm)
        │           ├── y (mm)
        │           └── z (mm)
        ├── General
        │   ├── authors
        │   ├── date
        │   ├── doi
        │   ├── original_filename
        │   ├── notes
        │   ├── time
        │   ├── time_zone
        │   └── title
        ├── Sample
        │   ├── credits
        │   ├── description
        │   ├── elements
        │   ├── thickness
        │   └── xray_lines
        └── Signal
            ├── FFT
            │   └── shifted
            ├── Noise_properties
            │   ├── Variance_linear_model
            │   │   ├── correlation_factor
            │   │   ├── gain_factor
            │   │   ├── gain_offset
            │   │   └── parameters_estimation_method
            │   └── variance
            ├── binned
            ├── quantity
            ├── signal_type
            └── signal_origin
    =#
    h5open(filename,"r") do hdf5
        ff = read_attributes(hdf5)
        @assert get(ff, "file_format", "") == "HyperSpy" "This file does not appear to be a HyperSpy file."
        ver = parse(Float64, get(ff, "file_format_version", 0.0))
        (ver < 1.2) && @warn "This reader may not be compatible with version $ver of HyperSpy data files."
        props = Dict{Symbol,Any}()
        props[:Filename] = filename
        if haskey(hdf5, "Experiments")
            exps = keys(hdf5["Experiments"])
            exp = something(name, first(exps))
            exp_key = "Experiments/$exp" 
            sig_key = "$exp_key/metadata/Signal"
            if haskey(hdf5, sig_key)
            #=  └── Signal
                    ├── FFT
                    │   └── shifted # N/A
                    ├── Noise_properties # N/A
                    │   ├── Variance_linear_model
                    │   │   ├── correlation_factor
                    │   │   ├── gain_factor
                    │   │   ├── gain_offset
                    │   │   └── parameters_estimation_method
                    │   └── variance
                    ├── binned
                    ├── quantity
                    ├── signal_type
                    └── signal_origin =#
                sig = read_attributes(hdf5[sig_key])
                st = get(sig, "signal_type", "Unknown")
                if st == "EDS_TEM"
                    props[:HSSignalType] = "EDS_TEM"
                elseif st == "EDS_SEM"
                    props[:HSSignalType] = "EDS_SEM"
                elseif st == "EELS"
                    props[:HSSignalType] = "EELS"
                else
                    @warn "Unsupported signal type $st which loading a HyperSpy file."
                    assign(props, :HSSignalType, st)
                end
                @assert get(sig, "recorded_by", "spectrum") == "spectrum" "Orderings other than `spectrum` are not supported"
                assign(props, :HSSignalOrigin, get(sig, "signal_origin", missing))
                assign(props, :HSQuantity, get(sig, "quantity",missing))
                assign(props, :HSBinned, get(sig, "binned",missing))
            end            
            sem_key = "$exp_key/metadata/Acquisition_instrument/SEM"
            if haskey(hdf5, sem_key)
                @assert get(props, :HSSignalType, "EDS_SEM")=="EDS_SEM"
            #=  ├── Acquisition_instrument
                │   ├── SEM
                │   │   ├── Detector
                │   │   │   ├── detector_type
                │   │   │   └── EDS
                │   │   │       ├── azimuth_angle (º)
                │   │   │       ├── elevation_angle (º)
                │   │   │       ├── energy_resolution_MnKa (eV)
                │   │   │       ├── live_time (s)
                │   │   │       └── real_time (s)
                │   │   ├── beam_current (nA)
                │   │   ├── beam_energy (keV)
                │   │   ├── probe_area (nm²)
                │   │   ├── convergence_angle (mrad)
                │   │   ├── magnification
                │   │   ├── microscope
                │   │   ├── Stage
                │   │   │   ├── rotation (º)
                │   │   │   ├── tilt_alpha (º)
                │   │   │   ├── tilt_beta (º)
                │   │   │   ├── x (mm)
                │   │   │   ├── y (mm)
                │   │   │   └── z (mm)
                │   │   └── working_distance (mm) =#
                assign(props, :InstrumentType, "SEM")
                sem_attr=read_attributes(hdf5[sem_key])
                assign(props, :ProbeCurrent, get(sem_attr, "beam_current", missing)) # nA
                assign(props, :BeamEnergy, get(sem_attr, "beam_energy", missing)*1000.0) # keV
                assign(props, :HSProbeArea, get(sem_attr, "probe_area", missing))
                assign(props, :HSConvergenceAngle, get(sem_attr, "convergence_angle", missing) * 0.001) # mrad
                assign(props, :HSMagnification, get(sem_attr, "magnification", missing))
                assign(props, :Instrument,  get(sem_attr, "microscope", missing))
                assign(props, :WorkingDistance, get(sem_attr,"working_distance", missing)*0.1) # mm
                if haskey(hdf5[sem_key],"Stage")
                    props[:Stage] = read_stage(read_attributes(hdf5["$sem_key/Stage"]))
                end
                det_key = "$sem_key/Detector"
                if haskey(hdf5, det_key)
                    det_attr = read_attributes(hdf5[det_key])
                    assign(props, :HSDetectorType, get(det_attr,"detector_type",missing))
                    eds_key = "$det_key/EDS"
                    if haskey(hdf5, eds_key)
                        read_eds(props, read_attributes(hdf5[eds_key]))
                    end
                end
            end
            tem_key = "$exp_key/metadata/Acquisition_instrument/TEM"
            if haskey(hdf5, tem_key)
            #=  ├── Acquisition_instrument
                │   └── TEM
                │       ├── Detector
                │       │   ├── EDS
                │       │   │   ├── azimuth_angle (º)
                │       │   │   ├── elevation_angle (º)
                │       │   │   ├── energy_resolution_MnKa (eV)
                │       │   │   ├── live_time (s)
                │       │   │   └── real_time (s)
                │       │   └── EELS
                │       │       ├── aperture (mm)
                │       │       ├── collection_angle (mrad)
                │       │       ├── dwell_time (s)
                │       │       ├── exposure (s)
                │       │       ├── frame_number
                │       │       └── spectrometer
                │       ├── Biprism
                │       │   ├── azimuth_angle (º)
                │       │   ├── position
                │       │   └── voltage (V)
                │       ├── acquisition_mode
                │       ├── beam_current (nA)
                │       ├── beam_energy (keV)
                │       ├── probe_area (nm²)
                │       ├── camera_length (mm)
                │       ├── convergence_angle (mrad)
                │       ├── magnification
                │       ├── microscope
                │       └── Stage
                │           ├── rotation (º)
                │           ├── tilt_alpha (º)
                │           ├── tilt_beta (º)
                │           ├── x (mm)
                │           ├── y (mm)
                │           └── z (mm) =#                
                assign(props, :InstrumentType, "TEM")
                tem_attr=read_attributes(hdf5[tem_key])
                assign(props, :ProbeCurrent, get(tem_attr, "beam_current", missing)) # nA
                assign(props, :BeamEnergy, get(tem_attr, "beam_energy", missing)*1000.0) # keV
                assign(props, :HSProbeArea, get(tem_attr, "probe_area", missing))
                assign(props, :HSConvergenceAngle, get(tem_attr, "convergence_angle", missing) * 0.001 )# mrad
                assign(props, :HSMagnification, get(tem_attr, "magnification", missing))
                assign(props, :Instrument,  get(tem_attr, "microscope", missing))
                assign(props, :HSCameraLength, get(tem_attr, "camera_length", missing)*0.1) # mm
                assign(props, :HSTEMTiltStage, get(tem_attr, "tilt_stage", missing) * deg2rad(1.0))
                assign(props, :HSAcquisitionMode, get(tem_attr, "acquisition_mode", missing))
                if haskey(hdf5[tem_key],"Stage")
                    props[:Stage] = read_stage(read_attributes(hdf5["$tem_key/Stage"]))
                end
                det_key = "$tem_key/Detector"
                if haskey(hdf5, det_key)
                    det_attr = read_attributes(hdf5[det_key])
                    assign(props, :HSDetectorType, get(det_attr,"detector_type",missing))
                    eds_key = "$det_key/EDS"
                    @assert get(props, :HSSignalType, "EDS_TEM")=="EDS_TEM"
                    if haskey(hdf5, eds_key)
                        read_eds(props, read_attributes(hdf5[eds_key]))
                    end
                    eels_key = "$det_key/EELS"
                    if haskey(hdf5, eels_key)
                        @assert get(props, :HSSignalType, "EELS")=="EELS"
                        eels_attr = read_attributes(hdf5[eels_attr])
                        assign(props, :HSAperture, get(eels_attr,"aperture",missing))
                        assign(props, :HSCollectionAngle, get(eels_attr,"collection_angle",missing) * 0.001)
                        assign(props, :HSDwellTime, get(eels_attr,"dwell_time",missing))
                        assign(props, :HSExposure, get(eels_attr,"exposure",missing))
                        assign(props, :HSFrameNumber, get(eels_attr,"frame_number",missing))
                        assign(props, :HSSpectrometer, get(eels_attr,"spectrometer",missing))
                    end
                end
            end
            #= │       ├── Biprism
               │       │   ├── azimuth_angle (º)
               │       │   ├── position
               │       │   └── voltage (V) =#
            bi_key = "$tem_key/Biprism"
            if haskey(hdf5, bi_key)
                bi_attr = read_attributes(hdf5[bi_key])
                assign(props, :HSBiprismAzAngle, get(bi_key, "azimuth_angle", missing) * deg2rad(1.0))
                assign(props, :HSBiprismPosition, get(bi_key, "position", missing))
                assign(props, :HSBiprismVoltage, get(bi_key, "voltage", missing))
            end
            #=  ├── General
                    │   ├── authors
                    │   ├── date
                    │   ├── doi
                    │   ├── original_filename
                    │   ├── notes
                    │   ├── time
                    │   ├── time_zone
                    │   └── title =#
            gen_key = "$exp_key/metadata/General"
            if haskey(hdf5, gen_key)
                general = read_attributes(hdf5[gen_key])
                assign(props, :Owner, get(general, "authors", missing))
                assign(props, :AcquisitionTime, read_datetime(general)) 
                assign(props, :Description, get(general, "notes", missing))
                assign(props, :Name, get(general, "title", missing))
                assign(props, :DOI, get(general, "doi", missing))
                assign(props, :Filename, get(general, "original_filename", missing))
                assign(props, :HSNotes, get(general, "notes", missing))
            end
            samp_key = "$exp_key/metadata/Sample"
            if haskey(hdf5, samp_key)
                #=  ├── Sample
                    │   ├── credits
                    │   ├── description
                    │   ├── elements
                    │   ├── thickness
                    │   └── xray_lines =#
                samp_attr = read_attributes(hdf5[samp_key])
                assign(props, :Sample, get(samp_attr, "description", missing))
                assign(props, :HSCredits, get(samp_attr, "credits", missing))
                elms = get(samp_attr, "elements", missing)
                if !ismissing(elms)
                    props[:Elements] = map(e->parse(Element, e), elms)
                end
                assign(props, :HSThickness, get(samp_attr, "thickness", missing))
                assign(props, :HSXRayLines, get(samp_attr, "xray_lines", missing))
            end
            scale, offset = 10.0, 0.0
            axisnames, fov, offsets = String[], Float64[], Float64[]
            i=0
            while haskey(hdf5["$exp_key"],"axis-$i")
                axis = read_attributes(hdf5["$exp_key/axis-$i"])
                unitname = uppercase(get(axis, "units", "eV")) 
                if unitname == "KEV"
                    unit = 1000.0
                elseif unitname == "EV"
                    unit = 1.0
                elseif unitname == "MM"
                    unit = 0.1
                elseif unitname == "NM"
                    unit = 1.0e-7
                else
                    unit = 1.0
                end
                push!(axisnames, uppercase(get(axis, "name", "axis-$i")))
                scale = get(axis, "scale", 10.0/unit) * unit
                offset = get(axis, "offset", 0.0) * unit
                push!(fov, scale * get(axis, "size", 1))
                push!(offsets, offset)
                i+=1
            end
            data = read(hdf5, "$exp_key/data")
            axisnames = reverse(axisnames[1:end-1])
            fov = reverse(fov[1:end-1])
            offsets = reverse(offsets[1:end-1])
            # Now swap the first two axes because column major vs row major...
            if length(axisnames)>1
                newdims = collect(1:length(axisnames)+1)
                newdims[2], newdims[3] = newdims[3], newdims[2]
                data = permutedims(data, newdims)
                axisnames[1], axisnames[2] = axisnames[2], axisnames[1]
                fov[1], fov[2] = fov[2], fov[1]
                offsets[1], offsets[2] = offsets[2], offsets[1]  
            end
            if haskey(props, :HSResolution)
                props[:Detector] = simpleEDS(size(data,1), scale, offset, props[:HSResolution])
            end
            return HyperSpectrum(LinearEnergyScale(offset, scale), props, data, #
                axisnames=axisnames, fov=fov, offset=offsets)
        end
    end
end

function write_hspy(filename::String, hspec::HyperSpectrum, mode=:SEM)
    h5open(filename) do hdf5
        write_hspy(hdf5, hspec, mode)
    end
end


function write_hspy(hdf5::File, hspec::HyperSpectrum, mode=:SEM)
    function assign_attr(gr, att, prop, tran=identity) 
        if haskey(hspec,prop)
            attributes(gr)[att] = tran(hspec[prop])
        end
    end
    groot = create_group(hdf5, "Experiments/$(hspec[:Name])")
    for (i, aa) in enumerate(AxisArrays.axes(h.counts))
        idx = ndims(h.counts)-i
        gax = create_group(groot, "axis-$idx")
        attributes(gax)["name"] = String(axisnames(h.counts)[idx])
        if i==1
            attributes(gax)["offset"] = first(aa.val)
            attributes(gax)["scale"] = step(aa.val)
            attributes(gax)["size"] = size(h.counts,i)
            attributes(gax)["units"] = "eV"
        else
            attributes(gax)["offset"] = first(aa.val)*1.0e7
            attributes(gax)["scale"] = step(aa.val)*1.0e7
            attributes(gax)["size"] = size(h.counts,i)
            attributes(gax)["units"] = "nm"
        end
    end
    geds = create_group(ginst, "Acquisition_instrument/"*String(mode)*"/Detector/EDS")
    assign_attr(geds, "azimuth_angle", :Azimuthal, rad2deg)
    assign_attr(geds, "elevation_angle", :TakeOffAngle, rad2deg)
    assign_attr(geds, "energy_resolution_MnKa", :Resolution)
    ggen = create_group(groot, "General")
    if haskey(attributes(gdet), :AcquisitionTime) 
        attributes(ggen)["date"] = "dd.mm.yyyy"
        attributes(ggen)["time"] = "hh.mm.ss"
        attributes(ggen)["time_zone"] = "+5.0"
    end
    gsamp = create_group(groot, "Sample")
    assign_attr(gsamp, "description", :Description)
    els = elms(hspec)
    if length(els)>0
        attributes(gsamp)["elements"] = join(symbol.(els),",")
    end
    gsamp = create_group(groot, "Signal")
    attributes(gsamp)["binned"] = true
    attributes(gsamp)["recorded_by"] = "spectrum"
    attributes(gsamp)["signal_type"] = "EDS_"*String(mode)
    groot["data"] = hspec.counts
end
    
"""
    ishspy(filename::String)::Bool

Is the file a HyperSpy-style HDF5 file?
"""
function ishspy(filename::String)::Bool
    return ishdf5(filename) && h5open(filename,"r") do h
        a = attributes(h)
        haskey(a,"file_format") && read(a["file_format"])=="HyperSpy"
    end
end