"""
Unfortunately, while the standard ImageMagick TIFF reader can read the image data from ASPEX TIFF files,
it can't read the spectrum data.  So I've written a lean TIFF tag reader to access the tags and
to read the spectral data into a Spectrum object.  This is then integrated into the FileIO mechanism
to permit using functions like `load(File(format"ASPEX TIFF"),filename)` to return a Spectrum
object loaded from the specified file.
"""

struct _ATField
    tagId::UInt16
    tagType::UInt16
    tagCount::UInt32
    tagData::Any

    function _ATField(ios, order)
        _TIFF_TYPES = (
            UInt8, # byte
            UInt8, # ASCII
            UInt16, # short
            UInt32, # long
            Rational{UInt32}, # rational
            Int8, # signed byte
            Int8, # undefined
            Int16, # signed short
            Int32, # signed long
            Rational{Int32}, # signed rational
            Float32, # float
            Float64, # double
        )
        tagId = read(ios, UInt16)
        tagType = read(ios, UInt16)
        tagCount = read(ios, UInt32)
        if tagCount * sizeof(_TIFF_TYPES[tagType]) <= 4
            tagData = [read(ios, _TIFF_TYPES[tagType]) for _ = 1:tagCount]
            extra = 4 - sizeof(_TIFF_TYPES[tagType]) * tagCount
            seek(ios, position(ios) + extra) # make sure that we have read 4 bytes
        else
            tagOffset = read(ios, UInt32)
            ret = position(ios)
            seek(ios, tagOffset)
            tagData = [read(ios, _TIFF_TYPES[tagType]) for _ = 1:tagCount]
            seek(ios, ret)
        end
        return new(tagId, tagType, tagCount, tagData)
    end
end

struct _TIFFIFD
    ifdTags::Vector{_ATField}
    ifdNext::UInt
    function _TIFFIFD(ios, order::Int)
        sz = read(ios, UInt16)
        ifdTags = [_ATField(ios, order) for _ = 1:sz]
        ifdNext = read(ios, UInt32)
        return new(ifdTags, ifdNext)
    end
end

function Base.get(ifd::_TIFFIFD, id::Integer, def)
    idx = findfirst(atf -> atf.tagId == id, ifd.ifdTags)
    return isnothing(idx) ? def : ifd.ifdTags[idx]
end

struct _TIFFInternals
    ifds::Vector{_TIFFIFD}
    function _TIFFInternals(ios)
        LITTLE_ENDIAN, BIG_ENDIAN, NOT_TIFF = 1, 2, 3
        magic = read(ios, UInt16)
        order = magic == 0x4D4D ? BIG_ENDIAN : (magic == 0x4949 ? LITTLE_ENDIAN : NOT_TIFF)
        fortytwo = read(ios, UInt16)
        if (order == NOT_TIFF) || (fortytwo ≠ 42)
            @error "This file is not a valid TIFF file."
        end
        if order ≠ LITTLE_ENDIAN
            @error "ASPEX Spectrum TIFF files are always little-endian ordered."
        end
        offset = read(ios, UInt32)
        ifds = _TIFFIFD[]
        while offset ≠ 0
            seek(ios, offset)
            ifd = _TIFFIFD(ios, order)
            push!(ifds, ifd)
            offset = ifd.ifdNext
        end
        return new(ifds)
    end
end

function _parseDesc(desc::AbstractString)
    number(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)", v)[1])
    date, time = missing, missing
    res, stagePos = Dict{Symbol,Any}(), Dict{Symbol,Float64}()
    for item in split(desc, "\n")
        (k, v) = split(item, "=")
        k = uppercase(k)
        try
            if k == "MAG"
                res[:ImageMag] = number(v)
                res[:FieldOfView] = (3.5 * 25.4 * 1.0e3) / number(v) # Assumed 3.5" image width in microns
            elseif k == "ZOOM"
                res[:ImageZoom] = number(v)
            elseif k == "ANALYSIS_DATE"
                date = Date(DateTime(v, "m/d/Y"))
            elseif k == "ANALYSIS_TIME"
                try
                    time = Time(DateTime(v, "H:M:S \\P\\M")) + Hour(12)
                catch
                    time = Time(DateTime(v, "H:M:S \\A\\M"))
                end # try
            elseif k == "OPERATOR"
                res[:Operator] = string(v)
            elseif k == "ACCELERATING_VOLTAGE"
                res[:BeamEnergy] = 1000.0 * number(v)
            elseif k == "WORKING_DISTANCE"
                res[:WorkingDistance] = 0.1 * number(v)
            elseif k == "COMMENT"
                res[:Name] = string(v)
            elseif k == "SAMPLE_DESCRIPTION"
                res[:Comment] = string(v)
            elseif k == "LIVE_TIME"
                res[:LiveTime] = number(v)
            elseif k == "ACQUISITION_TIME"
                res[:RealTime] = number(v)
            elseif k == "DETECTOR_TILT"
                res[:TakeOffAngle] = deg2rad(90.0 - number(v))
            elseif k == "STAGE_X"
                stagePos[:X] = 0.1 * number(v)
            elseif k == "STAGE_Y"
                stagePos[:Y] = 0.1 * number(v)
            elseif k == "STAGE_Z"
                stagePos[:Z] = 0.1 * number(v)
            elseif k == "STAGE_R"
                stagePos[:R] = deg2rad(number(v))
            elseif k == "STAGE_T"
                stagePos[:T] = deg2rad(number(v))
            elseif k == "STAGE_B"
                stagePos[:T] = deg2rad(number(v))
            end
        catch
            @info "Failed to parse $(k) => $(v)"
        end
    end
    if !isempty(stagePos)
        res[:StagePosition] = stagePos
    end
    if !ismissing(date)
        res[:AcquisitionTime] = ismissing(time) ? DateTime(date) : DateTime(date, time)
    end
    return res
end

function isAspexTIFF(filename)
    try
        open(filename, read = true) do ios
            return detectAspexTIFF(ios)
        end
    catch
        return false
    end
end

function detectAspexTIFF(ios)
    TIFF_SPECTRAL_DATA = 0x8352
    try
        seekstart(ios)
        ti = _TIFFInternals(ios)
        return any(ifd -> any(fld -> fld.tagId == TIFF_SPECTRAL_DATA, ifd.ifdTags), ti.ifds)
    catch
        return false
    end
end

function readAspexTIFF(
    file::AbstractString;
    withImgs = false,
    astype::Type{<:Real} = Float64,
)
    open(file, read = true) do ios
        return readAspexTIFF(ios, withImgs = withImgs, astype = astype)
    end
end

function readAspexTIFF(ios::IO; withImgs = false, astype::Type{<:Real} = Float64)
    # Special TIFF tags
    TIFF_SPECTRAL_DATA = 0x8352
    TIFF_SPECTRAL_XRES = 0x8353
    TIFF_SPECTRAL_XOFF = 0x8354
    TIFF_SPECTRAL_YRES = 0x8355
    TIFF_SPECTRAL_YOFF = 0x8356
    TIFF_IMAGE_DESCRIPTION = 270 # ASCII
    TIFF_SOFTWARE = 305 # ASCII
    floatonly(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)", v)[1])
    number(v) = parse(
        astype,
        match(
            astype isa Type{<:Integer} ? r"([+-]?[0-9]+)" : r"([+-]?[0-9]+[.]?[0-9]*)",
            v,
        )[1],
    )
    res = missing
    ti = _TIFFInternals(ios)
    for ifd in ti.ifds
        id = get(ifd, TIFF_IMAGE_DESCRIPTION, missing)
        sw = get(ifd, TIFF_SOFTWARE, missing)
        sp = get(ifd, TIFF_SPECTRAL_DATA, missing)
        sxr = get(ifd, TIFF_SPECTRAL_XRES, missing)
        sxo = get(ifd, TIFF_SPECTRAL_XOFF, missing)
        syr = get(ifd, TIFF_SPECTRAL_YRES, missing)
        syo = get(ifd, TIFF_SPECTRAL_YOFF, missing)
        if !ismissing(sp)
            @assert !(ismissing(sxr) || ismissing(sxo)) "X gain and offset data is missing from ASPEX TIFF spectrum"
            evperch = floatonly(ismissing(sxr) ? "10 eV/ch" : String(sxr.tagData))
            offset = floatonly(ismissing(sxo) ? "0 eV" : String(sxo.tagData))
            yoff = number(ismissing(syo) ? "0 counts" : String(syo.tagData))
            yres = number(ismissing(syr) ? "1 counts" : String(syr.tagData))
            energy = LinearEnergyScale(offset, evperch)
            data = map(i -> yoff + yres * convert(astype, i), sp.tagData)
            props = _parseDesc(String(id.tagData))
            if !ismissing(sw)
                props[:Instrument] = String(sw.tagData)
            end
            res = Spectrum(energy, data, props)
            break
        end
    end
    if withImgs && (!ismissing(res))
        try
            seekstart(ios)
            imgs = FileIO.load(Stream(format"TIFF", ios))
            nimgs = ndims(imgs)>2 ? size(imgs,3) : 1
            if haskey(res, :ImageMag) || haskey(res, :FieldOfView)
                fov = haskey(res, :ImageMag) ? (3.5*25.4)/res[:ImageMag] : 1.0e3*res[:FieldOfView] # X field-of-view in mm
                off = haskey(res, :StagePosition) ? (res[:StagePosition][:Y], res[:StagePosition][:X]) : (0.0, 0.0)
                ratio, pix = size(imgs,1) / size(imgs,2), fov / (size(imgs, 2)-1)
                ax = Axis{:x}((off[2]-0.5*fov)*mm:pix*mm:(off[2]-0.5*fov+pix*(size(imgs,2)-1))*mm)
                ay = Axis{:y}((off[1]+0.5*fov)*ratio*mm:-pix*mm:(off[1]+0.5*fov-pix*(size(imgs,1)-1))*ratio*mm)
                @assert length(ay)==size(imgs,1) && length(ax) == size(imgs,2)
                if nimgs == 2
                    # Macro image
                    res[Symbol("Image2")] = AxisArray(imgs[:, :, 2], ay, ax)
                    # Micro image
                    imgZoom = get(res, :ImageZoom, 1.0)  # >= 1.0
                    rfov, rpix = fov / imgZoom, fov / (imgZoom * (size(imgs, 2)-1))
                    ay = Axis{:y}(0.5*rfov*ratio*mm:-rpix*mm:(0.5*rfov-rpix*(size(imgs,1)-1))*ratio*mm)
                    ax = Axis{:x}(-0.5*rfov*mm:rpix*mm:(-0.5*rfov+rpix*(size(imgs,2)-1))*mm)
                    @assert length(ay)==size(imgs,1) && length(ax) == size(imgs,2)
                    res[Symbol("Image1")] = AxisArray(imgs[:, :, 1], ay, ax)
                else
                    # All images same FOV
                    if nimgs == 1
                        res[Symbol("Image1")] = AxisArray(imgs, ay, ax)
                    else
                        foreach(i -> res[Symbol("Image$i")] = AxisArray(imgs[:, :, i], ay, ax), 1:nimgs)
                    end
                end
            else
                # No scale data
                if nimgs == 1
                    res[Symbol("Image1")] = imgs
                else
                    foreach(i -> res[Symbol("Image$i")] = imgs[:, :, i], 1:nimgs)
                end
            end
        catch err
            @info err
            @info "Unable to read images from $(ios)."
        end
    end
    return res
end

"""
    LinearAlgebra.norm(aa::AxisArray, pt1, pt2, p::Real=2)

Returns the distance between `pt1` and `pt2` using the scaled axes in an AxisArray.
"""
LinearAlgebra.norm(aa::AxisArray, pt1, pt2, p::Real=2) = 
    norm( ( AxisArrays.axes(aa,i)[pt1[i]] for i in eachindex(pt1)) .- (AxisArrays.axes(aa,i)[pt2[i]] for i in eachindex(pt2)), p)
LinearAlgebra.norm(aa::AxisArray, pt1::CartesianIndex, pt2::CartesianIndex, p::Real=2) = 
    norm( ( AxisArrays.axes(aa,i)[pt1[i]] for i in Base.OneTo(length(pt1))) .- (AxisArrays.axes(aa,i)[pt2[i]] for i in Base.OneTo(length(pt1))), p)