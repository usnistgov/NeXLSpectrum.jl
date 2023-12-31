using TiffImages
using AxisArrays
import Unitful

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
    TIFF_SPECTRAL_DATA = UInt16(0x8352)
    try
        seekstart(ios)
        ti = TiffImages.load(ios; mmap=true)
        for ifd in ti.ifds
            if !ismissing(TiffImages.getdata(ifd, TIFF_SPECTRAL_DATA, missing))
                return true
            end
        end
        return false;
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
    TIFF_SPECTRAL_DATA = UInt16(0x8352)
    TIFF_SPECTRAL_XRES = UInt16(0x8353)
    TIFF_SPECTRAL_XOFF = UInt16(0x8354)
    TIFF_SPECTRAL_YRES = UInt16(0x8355)
    TIFF_SPECTRAL_YOFF = UInt16(0x8356)

    floatonly(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)", v)[1])
    number(v) = parse(
        astype,
        match(
            astype isa Type{<:Integer} ? r"([+-]?[0-9]+)" : r"([+-]?[0-9]+[.]?[0-9]*)",
            v,
        )[1],
    )
    res = missing
    ti = TiffImages.load(ios; mmap=true)
    for ifd in ti.ifds
        id = TiffImages.getdata(ifd, TiffImages.IMAGEDESCRIPTION, missing)
        sw = TiffImages.getdata(ifd, TiffImages.SOFTWARE, missing)
        sp = TiffImages.getdata(ifd, TIFF_SPECTRAL_DATA, missing)
        sxr = TiffImages.getdata(ifd, TIFF_SPECTRAL_XRES, missing)
        sxo = TiffImages.getdata(ifd, TIFF_SPECTRAL_XOFF, missing)
        syr = TiffImages.getdata(ifd, TIFF_SPECTRAL_YRES, missing)
        syo = TiffImages.getdata(ifd, TIFF_SPECTRAL_YOFF, missing)
        if !ismissing(sp)
            @assert !(ismissing(sxr) || ismissing(sxo)) "X gain and offset data is missing from ASPEX TIFF spectrum"
            evperch = floatonly(ismissing(sxr) ? "10 eV/ch" : String(sxr))
            offset = floatonly(ismissing(sxo) ? "0 eV" : String(sxo))
            yoff = number(ismissing(syo) ? "0 counts" : String(syo))
            yres = number(ismissing(syr) ? "1 counts" : String(syr))
            energy = LinearEnergyScale(offset, evperch)
            data = map(i -> yoff + yres * convert(astype, i), sp)
            props = _parseDesc(String(id))
            if !ismissing(sw)
                props[:Instrument] = String(sw)
            end
            res = Spectrum(energy, data, props)
        end
    end
    if withImgs && (!ismissing(res))
        try
            nimgs = ndims(ti)>2 ? size(ti,3) : 1
            if haskey(res, :ImageMag) || haskey(res, :FieldOfView)
                fov = haskey(res, :ImageMag) ? (3.5*25.4)/res[:ImageMag] : 1.0e3*res[:FieldOfView] # X field-of-view in mm
                off = haskey(res, :StagePosition) ? (res[:StagePosition][:Y], res[:StagePosition][:X]) : (0.0, 0.0)
                ratio, pix = size(ti,1) / size(ti,2), fov / (size(ti, 2)-1)
                ax = Axis{:x}((off[2]-0.5*fov)*Unitful.mm:pix*Unitful.mm:(off[2]-0.5*fov+pix*(size(ti,2)-1))*Unitful.mm)
                ay = Axis{:y}((off[1]+0.5*fov)*ratio*Unitful.mm:-pix*Unitful.mm:(off[1]+0.5*fov-pix*(size(ti,1)-1))*ratio*Unitful.mm)
                @assert length(ay)==size(ti,1) && length(ax) == size(ti,2)
                if nimgs == 2
                    # Macro image
                    res[Symbol("Image2")] = AxisArray(ti[:, :, 2], ay, ax)
                    # Micro image
                    imgZoom = get(res, :ImageZoom, 1.0)  # >= 1.0
                    rfov, rpix = fov / imgZoom, fov / (imgZoom * (size(ti, 2)-1))
                    ay = Axis{:y}(0.5*rfov*ratio*Unitful.mm:-rpix*Unitful.mm:(0.5*rfov-rpix*(size(ti,1)-1))*ratio*Unitful.mm)
                    ax = Axis{:x}(-0.5*rfov*Unitful.mm:rpix*Unitful.mm:(-0.5*rfov+rpix*(size(ti,2)-1))*Unitful.mm)
                    @assert length(ay)==size(ti,1) && length(ax) == size(ti,2)
                    res[Symbol("Image1")] = AxisArray(ti[:, :, 1], ay, ax)
                else
                    # All images same FOV
                    if nimgs == 1
                        res[Symbol("Image1")] = AxisArray(ti, ay, ax)
                    else
                        foreach(i -> res[Symbol("Image$i")] = AxisArray(ti[:, :, i], ay, ax), 1:nimgs)
                    end
                end
            else
                # No scale data
                if nimgs == 1
                    res[Symbol("Image1")] = ti
                else
                    foreach(i -> res[Symbol("Image$i")] = ti[:, :, i], 1:nimgs)
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