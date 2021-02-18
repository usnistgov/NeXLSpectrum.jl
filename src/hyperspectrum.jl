"""
   HyperSpectrum(arr::Array{T<:Real}, energy::EnergyScale, props::Array{Symbol, Any})

   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}; axisnames = ( "X", "Y", "Z", "A", "B", "C" ), fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ))
   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::AxisArray)
   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, dims::NTuple{<:Integer}, depth::Int, type::Type{Real}; axisnames = ( "X", "Y", "Z", "A", "B", "C" ), fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )

The HyperSpectrum struct represents a multi-dimensional array of Spectrum objects.  The dimension of a HyperSpectrum
may be 1 for a traverse or a line-scan, 2 for a spectrum image or higher for time-series of spectrum images or
multi-slice spectrum images.

The first constructor is used to create a HyperSpectrum from a raw Array of data.  The second to construct a HyperSpectrum
from another HyperSpectrum or an AxisArray.  The third from a description of the intended contents.

  * `axisnames`: A list of the names by which the axis can be referred
  * `fov`: The full width of the dimension in mm.

HyperSpectra differ from Array{Spectrum} in that the spectra in a HyperSpectrum must share properties like
:LiveTime and :ProbeCurrent.  HyperSpectrum objects can refer to line-scans (1D), spectrum images (2D), slice-and-view
(3D), time sequenced images (3D), or higher dimension spectrum images.

Internally, HyperSpectrum reinterpretes an Array{T<:Real, N+1} as an Array{Spectrum{T<:Real},N-1}.

HyperSpectrum objects can be read from a RPL/RAW file (using `readrplraw(filenamebase::AbstractString)`) but
can be constructed from any Array{<:Real}.

HyperSpectrum objects can be indexed using the standard Julia array idioms including a single integer index or a
CartesianIndex.  For example, to iterate over every spectrum in a HyperSpectrum

    % Construct a 20 × 20 spectrum image with 2048 channels of [0,255].
    hs = HyperSpectrum{es, props, (20,20), 2048, UInt8}
    for idx in eachindex(hs)
        spec = hs[idx]   % get a Spectrum representing the 2048 channels of data at idx
        spec[22] = 1     % Set the 22nd channel to 1
    end
"""
struct HyperSpectrum{T<:Real,N} <: AbstractArray{Spectrum{T},N}
    counts::AxisArray{T} # Organized as (ch, c, r) or (ch, x, y)
    index::CartesianIndices # Spectrum indices
    energy::EnergyScale
    properties::Dict{Symbol,Any}

    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}; axisnames = ( "X", "Y", "Z", "A", "B", "C" ), fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ))
        axes = [ Axis{:Channel}(1:size(arr,1)), #
                ( Axis{Symbol(axisnames[i-1])}(fov[i-1]*(-0.5:1.0/(size(arr,i)-1):0.5)) for i in 2:ndims(arr) )... ]
        new{eltype(arr),ndims(arr) - 1}(
            AxisArray(arr, axes...),
            CartesianIndices(size(arr)[2:end]),
            energy,
            props,
        )
    end
    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::AxisArray)
        new{eltype(arr),ndims(arr) - 1}(
            arr,
            CartesianIndices(size(arr)[2:end]),
            energy,
            props,
        )
    end
    function HyperSpectrum(
        energy::EnergyScale,
        props::Dict{Symbol,Any},
        dims::NTuple,
        depth::Int,
        type::Type{Real};
        axisnames = ( "X", "Y", "Z", "A", "B", "C" ), 
        fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
    )
        axes = [ Axis{:Channel}(1:depth), #
            ( Axis{Symbol(axisnames[i])}(fov[i]*(-0.5:1.0/(dims[i]-1):0.5)) for i in eachindex(dims) )... ]
        new{type,length(dims)}(
            AxisArray(zeros(type, ( depth, dims...)), axes),
            CartesianIndices(dims...),
            energy,
            props,
        )
    end
end

function Base.show(io::IO, m::MIME"text/plain", hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$(eltype(hss.counts)),$(ndims(hss))}[$(get(hss.properties,:Name,"Unnamed")), $(hss.energy), $(depth(hss)) ch]",
    )
end

function Base.show(io::IO, m::MIME"text/html", hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "<p>$sz HyperSpectrum{$(eltype(hss.counts)),$(ndims(hss))}[$(get(hss.properties,:Name,"Unnamed")), $(hss.energy), $(depth(hss)) ch]</p>",
    )
end

function Base.show(io::IO, hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$(eltype(hss.counts)),$(ndims(hss))}[$(get(hss.properties,:Name,"Unnamed")), $(hss.energy), $(depth(hss)) ch]",
    )
end

Base.eltype(::HyperSpectrum{T,N}) where {T<:Real,N} = Spectrum{T}
Base.length(hss::HyperSpectrum) = prod(size(hss))
Base.ndims(hss::HyperSpectrum) = ndims(hss.index)
Base.size(hss::HyperSpectrum) = size(hss.counts)[2:end]
Base.size(hss::HyperSpectrum, n) = size(hss.index, n)
Base.axes(hss::HyperSpectrum) = Base.axes(hss.index)
Base.axes(hss::HyperSpectrum, n) = Base.axes(hss.index, n)
Base.eachindex(hss::HyperSpectrum) = Base.OneTo(length(hss))
Base.stride(hss::HyperSpectrum, k::Int) = stride(hss.counts, k + 1)
Base.strides(hss::HyperSpectrum) = strides(hss.counts)[2:end]
depth(hss::HyperSpectrum) = size(hss.counts, 1)
name(hss::HyperSpectrum) = get(hss.properties, :Name, "HyperSpectrum")

function coordinate(hss::HyperSpectrum, idx::NTuple{N, <:Int}) where { N }
    av = axisvalues(hss.counts)
    pos = get(hss.properties, :StagePosition, Dict{Symbol,Float64}(:X=>0.0, :Y=>0.0))
    th = get(hss.properties, :ImageRotation, 0.0)
    rot = [ cos(th) -sin(th) ; sin(th) cos(th) ]
    xy = [ pos[:X], pos[:Y] ] + rot*[ av[i+1][ii] for (i,ii) in enumerate(idx) ]
    return Dict{Symbol,Float64}(:X=>xy[1], :Y=>xy[2])
end
coordinate(hss::HyperSpectrum, ci::CartesianIndex) = coordinate(hss, ci.I)

"""
    dose(hss::HyperSpectrum)

Returns the product of the probe current and the live-time.
"""
dose(hss::HyperSpectrum) = hss[:ProbeCurrent] * hss[:LiveTime]

function Base.getindex(hss::HyperSpectrum{T,N}, idx::Int)::Spectrum{T} where {T<:Real,N}
    props = copy(hss.properties)
    props[:Cartesian] = hss.index[idx]
    return Spectrum(hss.energy, counts(hss)[:, hss.index[idx]], props)
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    idx::Vararg{Int,N},
)::Spectrum{T} where {T<:Real,N}
    props = copy(hss.properties)
    props[:Cartesian] = CartesianIndex(idx...)
    props[:Name] = "$(name(hss))[$(props[:Cartesian])]"
    props[:StageCoorinate]=coordinate(hss, (idx...))
    return Spectrum(hss.energy, counts(hss)[:, idx...], props)
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    idx...
)::HyperSpectrum{T,N} where {T<:Real, N}
    props = copy(hss.properties)
    props[:Name] = "$(name(hss))[$idx]"
    return HyperSpectrum(hss.energy, props, hss.counts[:, idx...])
end

"""
    region(hss::HyperSpectrum{T, N}, ranges::AbstractRange...)::HyperSpectrum where {T<:Real,N}

Creates a view of a HyperSpectrum to represent the range of pixels within the `hss` HyperSpectrum.  Does
not copy the data or properties so any modifications to the region are also made to `hss`.
"""
function region(
    hss::HyperSpectrum{T,N},
    ranges::AbstractRange...,
)::HyperSpectrum{T,N} where {T<:Real,N}
    HyperSpectrum(hss.energy, hss.properties, counts(hss)[:, ranges...])
end

Base.getindex(hss::HyperSpectrum, sy::Symbol) = getindex(hss.properties, sy)
Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol) =
    setindex!(hss.properties, val, sy)

NeXLCore.energy(ch::Integer, hss::HyperSpectrum) = energy(ch, hss.energy)
channel(energy::Real, hss::HyperSpectrum) = channel(energy, hss.energy)
rangeofenergies(ch::Integer, hss::HyperSpectrum) =
    (energy(ch, hss.energy), energy(ch + 1, hss.energy))
properties(hss::HyperSpectrum)::Dict{Symbol,Any} = hss.properties

"""
    counts(hss::HyperSpectrum{T,N})::Array{T,N+1}

Creates type-friendly view of the counts data array.  Use of this function helps to avoid performance
penalties associated with boxing/unboxing the counts data.

    counts(hss::HyperSpectrum{T,N}, ci::CartesianIndex)::Vector{T}

Access the counts data associated with the pixel `ci`.

    counts(hss::HyperSpectrum{T,N}, ci::CartesianIndex, ch::Int)::T

Access the counts data at the pixel represented by `ci` and the channel represented by `ch`.
"""
counts(hss::HyperSpectrum{T,N}) where {T<:Real,N} = convert(Array{T,N + 1}, hss.counts)
counts(hss::HyperSpectrum, ci::CartesianIndex, ch::Int) = counts(hss)[ch, ci]
counts(hss::HyperSpectrum, ci::CartesianIndex) = counts(hss)[:, ci]

"""
    compress(hss::HyperSpectrum)

Returns a HyperSpectrum with smaller or equal storage space to `hss` without losing any infomation (except AbstractFloat
which compresses to Float32 with loss of precision).  Can change the storage type and/or reduce the depth of hss.
"""
function compress(hss::HyperSpectrum{T,N}) where {T<:Real,N}
    data = counts(hss)
    (minval, maxval) = extrema(data)
    maxi = mapreduce(
        ci -> something(findlast(c -> c != 0.0, data[:, ci]), 0),
        max,
        CartesianIndices(hss),
    )
    data = maxi < depth(hss) ? data[1:maxi, axes(hss.counts)[2:end]...] : data
    if eltype(T) in (Int16, Int32, Int64)
        for newtype in filter(ty -> sizeof(T) > sizeof(ty), (Int8, Int16, Int32))
            if minval >= typemin(newtype) && maxval <= typemax(newtype)
                return HyperSpectrum(hss.energy, copy(hss.properties), newtype.(data))
            end
        end
    elseif T in (UInt16, UInt32, UInt64)
        for newtype in filter(ty -> sizeof(T) > sizeof(ty), (UInt8, UInt16, UInt32))
            if maxval <= typemax(newtype)
                return HyperSpectrum(hss.energy, copy(hss.properties), newtype.(data))
            end
        end
    elseif T isa Float64
        if sizeof(T) < sizeof(Float32) &&
           minval >= typemin(Float32) &&
           maxval <= typemax(Float32)
            return HyperSpectrum(hss.energy, copy(hss.properties), Float32.(data))
        end
    end
    return maxi < depth(hss) ? HyperSpectrum(hss.energy, copy(hss.properties), data) : hss
end

"""
   plane(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the HyperSpectrum and is stored as a Float64 to ensure that not information is lost.
"""
function plane(
    hss::HyperSpectrum{T,N},
    chs::AbstractUnitRange{<:Integer},
    normalize = false,
) where {T<:Real,N}
    res = map(ci -> sum(Float64.(counts(hss)[chs, ci])), hss.index)
    if normalize
        res /= maximum(res)
    end
    return res
end

"""
   plane(hss::HyperSpectrum, ch::Int, normalize=false)

Extracts a single plane from a HyperSpectrum. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, ch::Int, normalize = false)
    res = counts(hss)[ch, hss.index]
    if normalize
        res /= maximum(res)
    end
    return res
end

"""
    maxpixel(hss::HyperSpectrum, filt=ci->true)
    maxpixel(hss::HyperSpectrum, cis::CartesianIndices, filt=ci->true)
    maxpixel(hss::HyperSpectrum, mask::BitArray)
    minpixel(hss::HyperSpectrum, filt=ci->true)
    minpixel(hss::HyperSpectrum, cis::CartesianIndices, filt=ci->true)
    minpixel(hss::HyperSpectrum, mask::BitArray)

Compute Bright's Max-Pixel derived signal for the entire HyperSpectrum or a rectanglar sub-region.
"""
maxpixel(hss::HyperSpectrum, filt::Function = ci -> true) =
    maxpixel(hss, CartesianIndices(hss), filt)
maxpixel(hss::HyperSpectrum, mask::BitArray) =
    maxpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function maxpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real,N}
    data, res = counts(hss), zeros(T, depth(hss))
    for ci in filter(ci -> filt(ci), cis), ch in Base.OneTo(size(res, 1))
        res[ch] = max(res[ch], data[ch, ci])
    end
    return Spectrum(hss.energy, res, copy(hss.properties))
end
minpixel(hss::HyperSpectrum, filt::Function = ci -> true) =
    minpixel(hss, CartesianIndices(hss), filt)
minpixel(hss::HyperSpectrum, mask::BitArray) =
    minpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function minpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real,N}
    data, res = counts(hss), fill(typemax(T), (depth(hss),))
    for ci in filter(ci -> filt(ci), cis), ch in Base.OneTo(size(res, 1))
        res[ch] = min(res[ch], data[ch, ci])
    end
    return Spectrum(hss.energy, res, copy(hss.properties))
end
avgpixel(hss::HyperSpectrum, filt::Function = ci -> true) =
    avgpixel(hss, CartesianIndices(hss), filt)
avgpixel(hss::HyperSpectrum, mask::BitArray) =
    avgpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function avgpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real,N}
    data, res, cx = counts(hss), zeros(Float64, depth(hss)), 0
    for ci in filter(ci -> filt(ci), cis)
        for ch in Base.OneTo(size(res, 1))
            res[ch] += data[ch, ci]
        end
        cx += 1
    end
    return Spectrum(hss.energy, res ./ max(1.0, cx), copy(hss.properties))
end

function sumcounts(hss::HyperSpectrum, cis::CartesianIndices = CartesianIndices(hss))
    data = counts(hss)
    return [sum(Int.(data[:, ci])) for ci in cis]
end

"""
    indexofmaxpixel(hss::HyperSpectrum, ch::Int) # at channel `ch`
    indexofmaxpixel(hss::HyperSpectrum) # all channels
    indexofmaxpixel(hss::HyperSpectrum, ch::Int, cis::CartesianIndices)
    indexofmaxpixel(hss::HyperSpectrum, cis::CartesianIndices)

Find the indices producing the maximum value in data[ch] or data[:] within 'cis' or full
spatial dimensions.
"""
function indexofmaxpixel(hss::HyperSpectrum, ch::Int, cis::CartesianIndices)
    data = counts(hss)
    maxidx, max = cis[1], data[ch, cis[1]]
    for idx in cis
        if data[ch, idx] > max
            maxidx = idx
            max = data[ch, idx]
        end
    end
    return maxidx
end

indexofmaxpixel(hss::HyperSpectrum, ch::Int) =
    indexofmaxpixel(hss, ch, CartesianIndices(hss))

indexofmaxpixel(hss::HyperSpectrum) = indexofmaxpixel(hss, CartesianIndices(hss))

function indexofmaxpixel(hss::HyperSpectrum{T,N}, cis::CartesianIndices) where {T<:Real,N}
    res, cix = zeros(T, depth(hss)), CartesianIndex[hss.index[1] for _ = 1:depth(hss)]
    data = counts(hss)
    for ci in cis, i in filter(i -> data[i, ci] > res[i], eachindex(res))
        res[i] = data[i, ci]
        cix[i] = ci
    end
    return cix
end

"""
    Base.sum(hss::HyperSpectrum{T, N}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, mask::Union{BitArray, Array{Bool}}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, filt::Function) where {T<:Real, N}

Compute a sum spectrum for all or a subset of the pixels in `hss`.
"""

function Base.sum(
    hss::HyperSpectrum{T,N},
    mask::Union{BitArray,Array{Bool}},
) where {T<:Real,N}
    @assert size(mask) == size(hss) "$(size(mask)) ≠ $(size(hss))"
    RT = T isa Int ? Int64 : Float64
    data = counts(hss)
    vals = mapreduce(
        ci -> data[:, ci],
        (a, b) -> a .+= b,
        filter(ci -> mask[ci], hss.index),
        init = zeros(RT, depth(hss)),
    )
    props = copy(hss.properties)
    if haskey(hss.properties, :LiveTime)
        props[:LiveTime] = count(mask) * hss[:LiveTime]
    end
    return Spectrum(hss.energy, vals, props)
end

function Base.sum(hss::HyperSpectrum{T,N}) where {T<:Real,N}
    RT = T isa Int ? Int64 : Float64
    data = counts(hss)
    vals = mapreduce(
        ci -> data[:, ci],
        (a, b) -> a .+= b,
        hss.index,
        init = zeros(RT, depth(hss)),
    )
    props = copy(hss.properties)
    if haskey(hss.properties, :LiveTime)
        props[:LiveTime] = length(hss) * hss[:LiveTime]
    end
    return Spectrum(hss.energy, vals, props)
end

function Base.sum(hss::HyperSpectrum{T,N}, filt::Function) where {T<:Real,N}
    RT = T isa Int ? Int64 : Float64
    data = counts(hss)
    vals = mapreduce(
        ci -> data[:, ci],
        (a, b) -> a .+= b,
        filter(ci -> filt(hss, ci), hss.index),
        init = zeros(RT, depth(hss)),
    )
    props = copy(hss.properties)
    if haskey(hss.properties, :LiveTime)
        props[:LiveTime] = count(ci -> filt(hss, ci), hss.index) * hss[:LiveTime]
    end
    return Spectrum(hss.energy, vals, props)
end

using Images

"""
    roiimages(hss::HyperSpectrum, achs::AbstractVector{<:AbstractUnitRange{<:Integer}})

Create an array of Gray images representing the intensity in each range of channels in
in `achs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
function roiimages(hss::HyperSpectrum, achs::AbstractVector{<:AbstractUnitRange{<:Integer}})
    ps = map(chs -> plane(hss, chs, false), achs)
    maxval = maximum(map(p -> maximum(p), ps))
    return map(p -> Gray.(convert.(N0f8, p / maxval)), ps)
end

"""
    roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer})

Create a count map from the specified contiguous range of channels.
"""
roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}) =
    Gray.(convert.(N0f8, plane(hss, chs, true)))

"""
    roiimage(hss::HyperSpectrum, cxr::CharXRay, n=5)

Create a count map for the specified characteristic X-ray.
"""
function roiimage(hss::HyperSpectrum, cxr::CharXRay, n = 5)
    center = channel(energy(cxr), hss.energy)
    return roiimage(hss, center-n:center+n)
end

"""
    roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, n=5)

Create an array of Gray images representing the intensity in each of the CharXRay lines
in `cxrs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
function roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, n = 5)
    achs = map(
        cxr -> channel(energy(cxr), hss.energy)-n:channel(energy(cxr), hss.energy)+n,
        cxrs,
    )
    return roiimages(hss, achs)
end
