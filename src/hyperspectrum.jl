"""
   HyperSpectrum(arr::Array{T<:Real}, energy::EnergyScale, props::Array{Symbol, Any})

A HyperSpectrum takes a Array{T<:Real, N} and converts it into an Array{Spectrum{T<:Real},N-1}.  In other words,
it makes it easy to interpet an hyperdimensional Array of numbers as an Array of spectra.
"""
struct HyperSpectrum{T<:Real,N} <: AbstractArray{Spectrum{T},N}
    counts::Array{T}
    index::CartesianIndices # Spectrum indices
    energy::EnergyScale
    properties::Dict{Symbol,Any}

    HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}) =
        new{eltype(arr),ndims(arr) - 1}(arr, CartesianIndices(size(arr)[2:end]), energy, props)

    HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, dims::NTuple, depth::Int, type::Type{Real}) =
        new{type,length(dims)}(zeros(type, dims...), CartesianIndices(dims...), energy, props)
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

Base.eltype(hss::HyperSpectrum{T,N}) where {T<:Real, N} = Spectrum{T}
Base.length(hss::HyperSpectrum) = prod(size(hss))
#Base.ndims(hss::HyperSpectrum) = ndims(hss.index)
Base.size(hss::HyperSpectrum) = size(hss.counts)[2:end]
Base.size(hss::HyperSpectrum, n) = size(hss.counts, n + 1)
Base.axes(hss::HyperSpectrum) = Base.axes(hss.index)
Base.axes(hss::HyperSpectrum, n) = Base.axes(hss.index, n)
Base.eachindex(hss::HyperSpectrum) = Base.OneTo(length(hss))
Base.stride(hss::HyperSpectrum, k::Int) = stride(hss.counts, k + 1)
Base.strides(hss::HyperSpectrum) = strides(hss.counts)[2:end]
depth(hss::HyperSpectrum) = size(hss.counts, 1)
dose(hss::HyperSpectrum) = hss[:ProbeCurrent] * hss[:LiveTime]

function Base.getindex(hss::HyperSpectrum{T,N}, idx::Int)::Spectrum{T} where {T<:Real,N}
    props = copy(hss.properties)
    props[:Cartesian] = hss.index[idx]
    return Spectrum(hss.energy, counts(hss)[:, hss.index[idx]], props)
end
function Base.getindex(hss::HyperSpectrum{T,N}, idx::Vararg{Int,N})::Spectrum{T} where {T<:Real,N}
    props = copy(hss.properties)
    props[:Cartesian] = CartesianIndex(idx...)
    return Spectrum(hss.energy, counts(hss)[:, idx...], props)
end
function region(hss::HyperSpectrum{T,N}, ranges::AbstractRange...)::HyperSpectrum{T,N} where {T<:Real,N}
    HyperSpectrum(hss.energy, hss.properties, counts(hss)[:, ranges...])
end
Base.getindex(hss::HyperSpectrum, sy::Symbol) = getindex(hss.properties, sy)
Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol) = setindex!(hss.properties, val, sy)

NeXLCore.energy(ch::Integer, hss::HyperSpectrum) = energy(ch, hss.energy)
channel(energy::Real, hss::HyperSpectrum) = channel(energy, hss.energy)
rangeofenergies(ch::Integer, hss::HyperSpectrum) = (energy(ch, hss.energy), energy(ch + 1, hss.energy))
properties(hss::HyperSpectrum)::Dict{Symbol,Any} = hss.properties
counts(hss::HyperSpectrum, ci::CartesianIndex, ch::Int) = counts(hss)[ch, ci]
counts(hss::HyperSpectrum{T,N}) where {T<:Real, N} = convert(Array{T,N+1}, hss.counts)

"""
    compress(hss::HyperSpectrum)

Returns a HyperSpectrum with smaller or equal storage space to `hss` without losing any infomation (except AbstractFloat
which compresses to Float32 with loss of precision).  Can change the storage type and/or reduce the depth of hss.
"""
function compress(hss::HyperSpectrum{T,N}) where {T<:Real, N}
    data = counts(hss)
    (minval, maxval) = extrema(data)
    maxi = mapreduce(ci -> findlast(c -> c != 0.0, data[:, ci]), max, CartesianIndices(hss))
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
        if sizeof(T) < sizeof(Float32) && minval >= typemin(Float32) && maxval <= typemax(Float32)
            return HyperSpectrum(hss.energy, copy(hss.properties), Float32.(data))
        end
    end
    return maxi < depth(hss) ? HyperSpectrum(hss.energy, copy(hss.properties), data) : hss
end

"""
   plane(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}, normalize = false)
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
    minpixel(hss::HyperSpectrum, filt=ci->true)
    minpixel(hss::HyperSpectrum, cis::CartesianIndices, filt=ci->true)

Compute Bright's Max-Pixel derived signal for the entire HyperSpectrum or a rectanglar sub-region.
"""
maxpixel(hss::HyperSpectrum, filt::Function=ci->true) = maxpixel(hss, CartesianIndices(hss), filt)
maxpixel(hss::HyperSpectrum, mask::BitArray) = maxpixel(hss, CartesianIndices(hss), ci->mask[ci])
function maxpixel(hss::HyperSpectrum{T,N}, cis::CartesianIndices, filt::Function=ci->true) where {T<:Real, N}
    data, res = counts(hss), zeros(T, depth(hss))
    for ci in filter(ci->filt(ci), cis), ch in Base.OneTo(size(res,1))
        res[ch] = max(res[ch], data[ch, ci])
    end
    return Spectrum(hss.energy, res, copy(hss.properties))
end
minpixel(hss::HyperSpectrum, filt::Function=ci->true) = minpixel(hss, CartesianIndices(hss), filt)
minpixel(hss::HyperSpectrum, mask::BitArray) = minpixel(hss, CartesianIndices(hss), ci->mask[ci])
function minpixel(hss::HyperSpectrum{T,N}, cis::CartesianIndices, filt::Function=ci->true) where {T<:Real, N}
    data, res = counts(hss), fill(typemax(T), (depth(hss), ))
    for ci in filter(ci->filt(ci), cis), ch in Base.OneTo(size(res,1))
        res[ch] = min(res[ch], data[ch, ci])
    end
    return Spectrum(hss.energy, res, copy(hss.properties))
end
avgpixel(hss::HyperSpectrum, filt::Function=ci->true) = avgpixel(hss, CartesianIndices(hss), filt)
avgpixel(hss::HyperSpectrum, mask::BitArray) = avgpixel(hss, CartesianIndices(hss), ci->mask[ci])
function avgpixel(hss::HyperSpectrum{T,N}, cis::CartesianIndices, filt::Function=ci->true) where {T<:Real, N}
    data, res, cx = counts(hss), zeros(Float64, depth(hss)), 0
    for ci in filter(ci->filt(ci), cis)
        for ch in Base.OneTo(size(res,1))
            res[ch] += data[ch, ci]
        end
        cx+=1
    end
    return Spectrum(hss.energy, res ./ max(1.0, cx), copy(hss.properties))
end

function sumcounts(hss::HyperSpectrum, cis::CartesianIndices=CartesianIndices(hss))
    data = counts(hss)
    return [ sum(Int.(data[:,ci])) for ci in cis ]
end

"""
    indexofmaxpixel(hss::HyperSpectrum, ch::Int) # at channel `ch`
    indexofmaxpixel(hss::HyperSpectrum) # all channels
    indexofmaxpixel(hss::HyperSpectrum, ch::Int, cis::CartesianIndices)
    indexofmaxpixel(hss::HyperSpectrum, cis::CartesianIndices)

Find the coordinates producing the maximum value in data[ch] or data[:] within 'cis' or full
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

indexofmaxpixel(hss::HyperSpectrum, ch::Int) = indexofmaxpixel(hss, ch, CartesianIndices(hss))

indexofmaxpixel(hss::HyperSpectrum) = indexofmaxpixel(hss, CartesianIndices(hss))

function indexofmaxpixel(hss::HyperSpectrum{T,N}, cis::CartesianIndices) where {T<:Real, N}
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

function Base.sum(hss::HyperSpectrum{T,N}, mask::Union{BitArray,Array{Bool}}) where {T<:Real,N}
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
    vals = mapreduce(ci -> data[:, ci], (a, b) -> a .+= b, hss.index, init = zeros(RT, depth(hss)))
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
roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}) = Gray.(convert.(N0f8, plane(hss, chs, true)))

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
    achs = map(cxr -> channel(energy(cxr), hss.energy)-n:channel(energy(cxr), hss.energy)+n, cxrs)
    return roiimages(hss, achs)
end
