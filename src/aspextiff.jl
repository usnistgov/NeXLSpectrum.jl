"""
Unfortunately, while the standard ImageMagick TIFF reader can read the image data from ASPEX TIFF files,
it can't read the spectrum data.  So I've written a lean TIFF tag reader to access the tags and
to read the spectral data into a Spectrum object.  This is then integrated into the FileIO mechanism
to permit using functions like `load(File(format"ASPEX TIFF"),filename)` to return a Spectrum
object loaded from the specified file.
"""

using NeXLSpectrum
using Dates
using FileIO
using Images
using ImageMagick

const _TIFF_TYPES = (
    ("byte", UInt8),
    ("ASCII", UInt8),
    ("short", UInt16),
    ("long", UInt32),
    ("rational", Rational{UInt32}),
    ("signed byte", Int8),
    ("undefined", Int8),
    ("signed short", Int16),
    ("signed long", Int32),
    ("signed rational", Rational{Int32}),
    ("float", Float32),
    ("doble", Float64),
)

TIFFType = Union{UInt8,UInt16,UInt32,Rational{UInt32},Int8,Int16,Int32,Rational{Int32},Float32,Float64}

const LITTLE_ENDIAN = 1
const BIG_ENDIAN = 2
const NOT_TIFF = 3

struct _ATField
    tagId::UInt16
    tagType::UInt16
    tagCount::UInt32
    tagData

    function _ATField(ios, order)
        tagId = read(ios, UInt16)
        tagType = read(ios, UInt16)
        tagCount = read(ios, UInt32)
        if tagCount * sizeof(_TIFF_TYPES[tagType][2]) <= 4
            tagData = [read(ios, _TIFF_TYPES[tagType][2]) for _ = 1:tagCount]
            extra = 4 - sizeof(_TIFF_TYPES[tagType][2]) * tagCount
            seek(ios, position(ios) + extra) # make sure that we have read 4 bytes
        else
            tagOffset = read(ios, UInt32)
            ret = position(ios)
            seek(ios, tagOffset)
            tagData = [read(ios, _TIFF_TYPES[tagType][2]) for _ = 1:tagCount]
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

function Base.get(ifd::_TIFFIFD, id, def)
    ff = findfirst(fld -> fld.tagId == convert(UInt16, id), ifd.ifdTags)
    return isnothing(ff) ? missing : ifd.ifdTags[ff]
end

struct _TIFFInternals
    ifds::Vector{_TIFFIFD}
    function _TIFFInternals(ios)
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

const TIFF_SPECTRAL_DATA = 0x8352
const TIFF_SPECTRAL_XRES = 0x8353
const TIFF_SPECTRAL_XOFF = 0x8354
const TIFF_SPECTRAL_YRES = 0x8355
const TIFF_SPECTRAL_YOFF = 0x8356
const TIFF_IMAGE_DESCRIPTION = 270 # ASCII
const TIFF_SOFTWARE = 305 # ASCII

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
    try
        seekstart(ios)
        ti = _TIFFInternals(ios)
        return any(ifd->any(fld->fld.tagId==TIFF_SPECTRAL_DATA, ifd.ifdTags), ti.ifds)
    catch
        return false
    end
end

function readAspexTIFF(file::AbstractString; withImgs = false, astype::Type{<:Real} = Float64)
    open(file, read = true) do ios
        return readAspexTIFF(ios, withImgs = withImgs, astype = astype)
    end
end

function readAspexTIFF(ios::IO; withImgs = false, astype::Type{<:Real} = Float64)
    floatonly(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)", v)[1])
    number(v) = parse(astype, match(astype isa Type{<:Integer} ? r"([+-]?[0-9]+)" : r"([+-]?[0-9]+[.]?[0-9]*)", v)[1])
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
            res[:Image] = FileIO.load(Stream(format"TIFF", ios))
        catch err
            @info "Unable to read images from $(ios)."
        end
    end
    return res
end

const ASPEX_TIFF = format"ASPEX TIFF"

load(ios::Stream{ASPEX_TIFF}; withImgs = false) = readAspexTIFF(ios, withImgs)

function load(file::File{ASPEX_TIFF}; withImgs = false)
    open(filename) do ios
        return readAspexTIFF(ios)
    end
end


function save(f::Stream{ASPEX_TIFF}, data)
    @error "Saving to ASPEX TIFF streams is not implemented. Probably never will be."
end

function save(f::File{ASPEX_TIFF}, data)
    @error "Saving to ASPEX TIFF files is not implemented. Probably never will be."
end

FileIO.add_format(ASPEX_TIFF, detectAspexTIFF, [".tif", ".tiff"])
