using EzXML

"""
    readbrukerspx(fn::AbstractString)::Spectrum
    readbrukerspx(io::IO)::Spectrum

Read a Bruker SPX EDS spectrum file into a `Spectrum`.
"""
function readbrukerspx(filename::String)::Spectrum
    return open(filename,read=true) do io
        readbrukerspx(io)
    end
end

function readbrukerspx(io::IO)::Spectrum
    xml = readxml(io)
    # Read the counts data
    chdata = findall("//TRTSpectrum/ClassInstance/Channels",xml)[1].content
    counts = parse.(Int, strip.(split(chdata,",")))
    # Read the header data
    props = Dict{Symbol, Any}()
    nrgy = missing
    item = findfirst("//TRTSpectrum/ClassInstance/ClassInstance[@Type='TRTSpectrumHeader']",xml)
    if !isnothing(item)
        dmy = Date(DateTime(findfirst("Date",item).content,"d.m.Y")) # <Date>10.5.2017</Date>
        hms = Time(DateTime(findfirst("Time",item).content,"H:M:S")) # <Time>11:36:22</Time>
        props[:AcquisitionTime] = DateTime(dmy, hms)
        @assert parse(Int, findfirst("ChannelCount",item).content)==length(counts)
        off = 1000.0*parse(Float64,findfirst("CalibAbs",item).content) # <CalibAbs>-9.5639563E-1</CalibAbs>
        gain = 1000.0*parse(Float64,findfirst("CalibLin",item).content) # <CalibLin>1.0001E-2</CalibLin>
        nrgy = LinearEnergyScale(off, gain)
    end
    item = findfirst("//TRTSpectrum/ClassInstance/ClassInstance[@Type='TRTPSEElementList']",xml)
    if !isnothing(item)
        data = findall("ChildClassInstances/ClassInstance",item)
        elmv = Element[PeriodicTable.elements[parse(Int, findfirst("Element",datum).content)] for datum in data]
        if length(elmv)>0
            props[:Elements] = elmv
        end
    end
    # item = findall("//TRTSpectrum/ClassInstance/ClassInstance[@Type='TRTResult']",xml) # Someday...
    item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTSpectrumHardwareHeader']",xml)
    item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTSpectrumHardwareHeader']",xml)
    if !isnothing(item)
        props[:RealTime]=0.001*parse(Int, findfirst("RealTime",item).content) # <RealTime>310610</RealTime>
        props[:LiveTime]=0.001*parse(Int, findfirst("LifeTime",item).content) # <LifeTime>300294</LifeTime>
        props[:BrukerThroughput] = parse(Int, findfirst("ShapingTime", item).content) # <ShapingTime>130000</ShapingTime>
    end
    item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTDetectorHeader']",xml)
    if !isnothing(item)
        props[:DetectorSerialNumber] = findfirst("Serial",item).content
        props[:DetectorModel] = "Bruker "*findfirst("Type",item).content
        thickness = 0.1 * parse(Float64, findfirst("DetectorThickness",item).content) # <DetectorThickness>0.45</DetectorThickness>
        props[:DetectorThickness] = thickness
        sideadlayer = parse(Float64, findfirst("SiDeadLayerThickness", item).content) # <SiDeadLayerThickness>0.029</SiDeadLayerThickness>
        props[:DeadLayerThickness] = 0.1 * sideadlayer
        wls = findfirst("WindowLayers", item)
        window = Film[]
        for i in 0:100
            wl = findfirst("Layer$i",wls)
            if isnothing(wl)
                break
            end
            z = parse(Int, wl["Atom"])  #   <Layer0 Atom="4" Thickness="1.3E1"/>
            thk = 1.0e-7 * parse(Float64, wl["Thickness"]) # nm
            relarea = haskey(wl,"RelativeArea") ? parse(Float64, wl["RelativeArea"]) : 1.0
            if relarea == 1.0
                push!(window, Film(pure(PeriodicTable.elements[z]), thk))
            end
        end
        if length(window)>0
            props[:Window] = window
        end
    end
    item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTXrfHeader']",xml)
    if !isnothing(item)
        props[:BeamEnergy] = 1000.0 * parse(Float64,findfirst("Voltage",item).content) # <Voltage>50</Voltage>
        props[:XRFTubeAnode] = PeriodicTable.elements[parse(Int,findfirst("Anode",item).content)]  # <Anode>42</Anode>
        props[:ProbeCurrent] = 1.0e3*parse(Float64, findfirst("Current", item).content) # <Current>99</Current> ???microamps???
        props[:XRFTubeIncidentAngle] = deg2rad(parse(Float64, findfirst("TubeIncidentAngle", item).content)) # <TubeIncidentAngle>8.4E1</TubeIncidentAngle>
        props[:XRFTubeTakeOffAngle] = deg2rad(parse(Float64, findfirst("TubeTakeOffAngle", item).content)) # <TubeTakeOffAngle>6</TubeTakeOffAngle>
        props[:XRFExcitationAngle] = deg2rad(parse(Float64, findfirst("ExcitationAngle", item).content)) # <ExcitationAngle>5E1</ExcitationAngle>
        props[:Elevation] = deg2rad(parse(Float64, findfirst("DetectionAngle", item).content)) # <DetectionAngle>5E1</DetectionAngle>
        props[:XRFExcitationPathLength] = 0.1 * parse(Float64, findfirst("ExcitationPathLength", item).content) # <ExcitationPathLength>1E1</ExcitationPathLength>
        props[:XRFDetectionPathLength] = 0.1 * parse(Float64, findfirst("DetectionPathLength", item).content) # <DetectionPathLength>2E1</DetectionPathLength>
        props[:DetectorSolidAngle] = parse(Float64, findfirst("SolidAngleDetection", item).content)
        # <AzimutAngleAbs>0</AzimutAngleAbs>
        # <DetAzimutAngle>0</DetAzimutAngle>
        props[:ChamberPressure] = parse(Float64, findfirst("ChamberPressure", item).content) # <ChamberPressure>2E1</ChamberPressure>
        props[:ChamberAtmosphere] = findfirst("Atmosphere", item).content # <Atmosphere>Air</Atmosphere>
        props[:XRFSampleTilt] = deg2rad(parse(Float64, findfirst("TiltAngle", item).content)) # <TiltAngle>0</TiltAngle>
        props[:TakeOffAngle] = props[:Elevation] + props[:XRFSampleTilt]
        welm = PeriodicTable.elements[parse(Int,findfirst("TubeWindow/AtomicNumber",item).content)]
        wthk = 1.0e-4 * parse(Float64,findfirst("TubeWindow/Thickness",item).content)
        props[:XRFTubeWindow] = Film(pure(welm),wthk)
    end
    item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTESMAHeader']",xml)
    if !isnothing(item)
        props[:BeamEnergy] = 1000.0 * parse(Float64, findfirst("PrimaryEnergy", item).content)  # <PrimaryEnergy>2E1</PrimaryEnergy>
        props[:Elevation] = (props[:TakeOffAngle] = deg2rad(parse(Float64, findfirst("ElevationAngle", item).content)))
        props[:WorkingDistance] = 0.1 * parse(Float64, findfirst("WorkingDistance", item).content)
    end
    #item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTXrfFPModelHeader']",xml)
    #item = findfirst("//TRTSpectrum/ClassInstance/TRTHeaderedClass/ClassInstance[@Type='TRTQuantitatorConfig']",xml)
    return Spectrum(nrgy, counts, props)
end


function detectbrukerspx(io::IO)
    try
        startswith(strip(readline(io)),"<?xml version=\"1.0\" encoding=") &&
            strip(readline(io))=="<TRTSpectrum>";
    catch
        return false
    end
end
