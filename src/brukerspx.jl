# Read the Bruker XML-based EDS spectrum file format.

"""
    readbrukerspx(fn::AbstractString)::Spectrum
    readbrukerspx(io::IO)::Spectrum

Read a Bruker SPX EDS spectrum file into a `Spectrum`.
"""
function readbrukerspx(filename::String)::Spectrum
    return open(filename, read = true) do io
        readbrukerspx(io)
    end
end

function readbrukerspx(io::IO)::Spectrum
    rootxml = readxml(io)
    spxml = findfirst("//ClassInstance[@Type='TRTSpectrum']", rootxml)
    brukerxml2spectrum(spxml)
end

function brukerxml2spectrum(xml::Node)::Spectrum
    # Read the counts data
    chdata = findfirst("//Channels", xml).content
    counts = parse.(Int, strip.(split(chdata, ",")))
    # Read the header data
    props = Dict{Symbol,Any}()
    nrgy = missing
    item = findfirst(
        "//ClassInstance[@Type='TRTSpectrumHeader']",
        xml,
    )
    if !isnothing(item)
        dmy = Date(DateTime(findfirst("Date", item).content, "d.m.Y")) # <Date>10.5.2017</Date>
        hms = Time(DateTime(findfirst("Time", item).content, "H:M:S")) # <Time>11:36:22</Time>
        props[:AcquisitionTime] = DateTime(dmy, hms)
        @assert parse(Int, findfirst("ChannelCount", item).content) == length(counts)
        off = 1000.0 * parse(Float64, findfirst("CalibAbs", item).content) # <CalibAbs>-9.5639563E-1</CalibAbs>
        gain = 1000.0 * parse(Float64, findfirst("CalibLin", item).content) # <CalibLin>1.0001E-2</CalibLin>
        nrgy = LinearEnergyScale(off, gain)
        props[:SigmaZerokV] = 1000.0 * parse(Float64, findfirst("SigmaAbs", item).content)  # useful for ZeroPeak Integration
    end
    item = findfirst(
        "//ClassInstance[@Type='TRTPSEElementList']",
        xml,
    )
    if !isnothing(item)
        data = findall("ChildClassInstances/ClassInstance", item)
        elmv = Element[
            PeriodicTable.elements[parse(Int, findfirst("Element", datum).content)] for
            datum in data
        ]
        if length(elmv) > 0
            props[:Elements] = elmv
        end
    end
    # item = findall("//TRTSpectrum/ClassInstance/ClassInstance[@Type='TRTResult']",xml) # Someday...
    item = findfirst(
        "//TRTHeaderedClass/ClassInstance[@Type='TRTSpectrumHardwareHeader']",
        xml,
    )
    if !isnothing(item)
        if !isnothing(findfirst("RealTime", item))
            props[:RealTime] = 0.001 * parse(Int, findfirst("RealTime", item).content) # <RealTime>310610</RealTime>
        end
        if !isnothing(findfirst("LifeTime", item))
            props[:LiveTime] = 0.001 * parse(Int, findfirst("LifeTime", item).content) # <LifeTime>300294</LifeTime>
        end
        props[:BrukerThroughput] = parse(Int, findfirst("ShapingTime", item).content) # <ShapingTime>130000</ShapingTime>
        # in case Real time and LifeTime is not there (Hyperspectra) it can
        # be calculated from ZeroPeak knowing following parameters: 
        props[:ZeroPeakPosition] = parse(Int, findfirst("ZeroPeakPosition", item).content)  # in Channel index (0 or 1 based? Delfi Pascal has both)
        props[:ZeroPeakFrequency] = parse(Int, findfirst("ZeroPeakFrequency", item).content)  # in Hz
    end
    item = findfirst(
        "//TRTHeaderedClass/ClassInstance[@Type='TRTDetectorHeader']",
        xml,
    )
    if !isnothing(item)
        props[:DetectorSerialNumber] = findfirst("Serial", item).content
        props[:DetectorModel] = "Bruker " * findfirst("Type", item).content
        thickness = 0.1 * parse(Float64, findfirst("DetectorThickness", item).content) # <DetectorThickness>0.45</DetectorThickness>
        props[:DetectorThickness] = thickness
        sideadlayer = parse(Float64, findfirst("SiDeadLayerThickness", item).content) # <SiDeadLayerThickness>0.029</SiDeadLayerThickness>
        props[:DeadLayerThickness] = 0.1 * sideadlayer
        wls = findfirst("WindowLayers", item)
        window = Film[]
        for i = 0:100
            wl = findfirst("Layer$i", wls)
            if isnothing(wl)
                break
            end
            z = parse(Int, wl["Atom"])  #   <Layer0 Atom="4" Thickness="1.3E1"/>
            thk = 1.0e-7 * parse(Float64, wl["Thickness"]) # nm
            relarea = haskey(wl, "RelativeArea") ? parse(Float64, wl["RelativeArea"]) : 1.0
            if relarea == 1.0
                push!(window, Film(pure(PeriodicTable.elements[z]), thk))
            end
        end
        if length(window) > 0
            props[:Window] = window
        end
    end
    item = findfirst(
        "//TRTHeaderedClass/ClassInstance[@Type='TRTXrfHeader']",
        xml,
    )
    if !isnothing(item)
        # <Voltage>50</Voltage>
        props[:BeamEnergy] = 1000.0 * parse(Float64, findfirst("Voltage", item).content)
        # <Anode>42</Anode> 
        props[:XRFTubeAnode] =
            PeriodicTable.elements[parse(Int, findfirst("Anode", item).content)]
        # <Current>99</Current> ???microamps???
        props[:ProbeCurrent] = 1.0e3 * parse(Float64, findfirst("Current", item).content)
        # <TubeIncidentAngle>8.4E1</TubeIncidentAngle> 
        props[:XRFTubeIncidentAngle] =
            deg2rad(parse(Float64, findfirst("TubeIncidentAngle", item).content))
        # <TubeTakeOffAngle>6</TubeTakeOffAngle>
        props[:XRFTubeTakeOffAngle] =
            deg2rad(parse(Float64, findfirst("TubeTakeOffAngle", item).content))
        # <ExcitationAngle>5E1</ExcitationAngle>
        props[:XRFExcitationAngle] =
            deg2rad(parse(Float64, findfirst("ExcitationAngle", item).content))
        # <DetectionAngle>5E1</DetectionAngle>
        props[:Elevation] =
            deg2rad(parse(Float64, findfirst("DetectionAngle", item).content))
        # <ExcitationPathLength>1E1</ExcitationPathLength>
        props[:XRFExcitationPathLength] =
            0.1 * parse(Float64, findfirst("ExcitationPathLength", item).content)
        # <DetectionPathLength>2E1</DetectionPathLength>
        props[:XRFDetectionPathLength] =
            0.1 * parse(Float64, findfirst("DetectionPathLength", item).content)
        props[:DetectorSolidAngle] =
            parse(Float64, findfirst("SolidAngleDetection", item).content)
        # <AzimutAngleAbs>0</AzimutAngleAbs>
        # <DetAzimutAngle>0</DetAzimutAngle>
        # <ChamberPressure>2E1</ChamberPressure>
        props[:ChamberPressure] = parse(Float64, findfirst("ChamberPressure", item).content)
        # <Atmosphere>Air</Atmosphere>
        props[:ChamberAtmosphere] = findfirst("Atmosphere", item).content
        # <TiltAngle>0</TiltAngle>
        props[:XRFSampleTilt] =
            deg2rad(parse(Float64, findfirst("TiltAngle", item).content))
        props[:TakeOffAngle] = props[:Elevation] + props[:XRFSampleTilt]
        welm = PeriodicTable.elements[parse(
            Int,
            findfirst("TubeWindow/AtomicNumber", item).content,
        )]
        wthk = 1.0e-4 * parse(Float64, findfirst("TubeWindow/Thickness", item).content)
        props[:XRFTubeWindow] = Film(pure(welm), wthk)
    end
    item = findfirst(
        "//TRTHeaderedClass/ClassInstance[@Type='TRTESMAHeader']",
        xml,
    )
    if !isnothing(item)
        # <PrimaryEnergy>2E1</PrimaryEnergy>
        props[:BeamEnergy] =
            1000.0 * parse(Float64, findfirst("PrimaryEnergy", item).content)
        props[:Elevation] = (
            props[:TakeOffAngle] =
                deg2rad(parse(Float64, findfirst("ElevationAngle", item).content))
        )
        wd = findfirst("WorkingDistance", item)
        if !isnothing(wd)
            props[:WorkingDistance] = 0.1 * parse(Float64, wd.content)
        end
    end
    #item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTXrfFPModelHeader']",xml)
    #item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTQuantitatorConfig']",xml)
    return Spectrum(nrgy, counts, props)
end


function detectbrukerspx(io::IO)
    try
        startswith(strip(readline(io)), "<?xml version=\"1.0\" encoding=") &&
            strip(readline(io)) == "<TRTSpectrum>"
    catch
        return false
    end
end
