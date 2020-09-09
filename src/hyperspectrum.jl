"""
   HyperSpectrum(arr::Array{T<:Real}, energy::EnergyScale, props::Array{Symbol, Any})

A HyperSpectrum takes a Array{T<:Real, N} and converts it into an Array{Spectrum{T<:Real},N-1}.  In other words,
it makes it easy to interpet an hyperdimensional Array of numbers as an Array of spectra.
"""
struct HyperSpectrum{T<:Real, N} <: AbstractArray{Spectrum{T}, N}
    counts::Array{T}
    index::CartesianIndices # Spectrum indices
    energy::EnergyScale
    properties::Dict{Symbol, Any}

    HyperSpectrum(energy::EnergyScale, props::Dict{Symbol, Any}, arr::Array{<:Real}) =
        new{eltype(arr), ndims(arr)-1}(arr, CartesianIndices(size(arr)[2:end]), energy, props)

    HyperSpectrum(energy::EnergyScale,
            props::Dict{Symbol, Any},
            dims::NTuple,
            depth::Int,
            type::Type{Real}
        ) =
        new{type, length(dims)}(zeros(type, dims...), CartesianIndices(dims...), energy, props)
end

function Base.show(io::IO, hss::HyperSpectrum{T,N}) where {T<:Real,N}
    print(io,"HyperSpectrum{$T,$N}[$(get(hss.properties,:Name,"Unnamed")), $(hss.energy), $(size(hss.index))]")
end


Base.eltype(hss::HyperSpectrum) = Spectrum{eltype(hss.counts)}
Base.length(hss::HyperSpectrum) = prod(size(hss))
#Base.ndims(hss::HyperSpectrum) = ndims(hss.index)
Base.size(hss::HyperSpectrum) = size(hss.counts)[2:end]
Base.size(hss::HyperSpectrum, n) = size(hss.counts, n+1)
Base.axes(hss::HyperSpectrum) = Base.axes(hss.index)
Base.axes(hss::HyperSpectrum, n) = Base.axes(hss.index, n)
Base.eachindex(hss::HyperSpectrum) = Base.OneTo(length(hss))
Base.stride(hss::HyperSpectrum, k::Int) = stride(hss.counts, k+1)
Base.strides(hss::HyperSpectrum) = strides(hss.counts)[2:end]
depth(hss::HyperSpectrum) = size(hss.counts,1)
dose(hss::HyperSpectrum) = hss[:ProbeCurrent]*hss[:LiveTime]

function Base.getindex(hss::HyperSpectrum{T,N}, idx::Int)::Spectrum{T} where {T<:Real, N}
    props = copy(hss.properties)
    props[:Cartesian] = hss.index[idx]
    return Spectrum(hss.energy, hss.counts[:, hss.index[idx]], props)
end
function Base.getindex(hss::HyperSpectrum{T,N}, idx::Vararg{Int, N})::Spectrum{T} where {T<:Real, N}
    props = copy(hss.properties)
    props[:Cartesian] = CartesianIndex(idx...)
    return Spectrum(hss.energy, hss.counts[:, idx...], props)
end
function region(hss::HyperSpectrum{T,N}, ranges::AbstractRange...)::HyperSpectrum{T,N} where {T<:Real, N}
    HyperSpectrum(hss.energy, hss.properties, hss.counts[:, ranges...])
end
Base.getindex(hss::HyperSpectrum, sy::Symbol) =
    getindex(hss.properties, sy)
Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol) =
    setindex!(hss.properties, val, sy)

NeXLCore.energy(ch::Integer, hss::HyperSpectrum) = energy(ch, hss.energy)
channel(energy::Real, hss::HyperSpectrum) = channel(energy, hss.energy)
rangeofenergies(ch::Integer, hss::HyperSpectrum) =
    ( energy(ch, hss.energy), energy(ch+1, hss.energy) )
properties(hss::HyperSpectrum)::Dict{Symbol, Any} = hss.properties
counts(hss::HyperSpectrum, ci::CartesianIndex, ch::Int) = hss.counts[ch,ci]
"""
    compressed(hss::HyperSpectrum)

Returns a HyperSpectrum with smaller or equal storage space to `hss` without losing any infomation (except Float64 which
compresses to Float32 with loss of precision).
"""
function compressed(hss::HyperSpectrum)
    maxval = maximum(hss.counts)
    if (eltype(hss.counts) in ( UInt16, UInt32 )) && (maxval <= 255)
        return HyperSpectrum(UInt8.(hss.counts), hss.energy, hss.properties)
    elseif (eltype(hss.counts) isa UInt32 ) && (maxval <= 65535)
        return HyperSpectrum(UInt16.(hss.counts), hss.energy, hss.properties)
    elseif (eltype(hss.counts) in ( Int16, Int32 )) && (maxval <= 127)
        return HyperSpectrum(Int8.(hss.counts), hss.energy, hss.properties)
    elseif (eltype(hss.counts) isa Int32 ) && (maxval <= 32767)
        return HyperSpectrum(Int16.(hss.counts), hss.energy, hss.properties)
    elseif eltype(hss.counts) isa Float64
        return HyperSpectrum(Float32.(hss.counts), hss.energy, hss.properties)
    end
    return hss
end

"""
   plane(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer}, normalize=false)
    res = map(ci->sum(Float64.(hss.counts[chs,ci])), hss.index)
    if normalize
        res/=maximum(res)
    end
    return res
end

"""
   plane(hss::HyperSpectrum, ch::Int, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, ch::Int, normalize=false)
    res = hss.counts[ch,hss.index]
    if normalize
        res/=maximum(res)
    end
    return res
end

"""
    maxpixel(hss::HyperSpectrum)

Compute Bright's Max-Pixel derived signal.
"""
function maxpixel(hss::HyperSpectrum)
    maxi(i) = maximum( ( hss.counts[idx] for idx in i:stride(hss.counts,2):length(hss.counts) ))
    return Spectrum(hss.energy, map(i->maxi(i), Base.axes(hss.counts,1)), hss.properties)
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
    maxidx, max = cis[1], hss.counts[ch, cis[1]]
    for idx in cis
        if hss.counts[ch, idx] > max
            maxidx = idx
            max=hss.counts[ch, idx]
        end
    end
    return maxidx
end

indexofmaxpixel(hss::HyperSpectrum, ch::Int) =
    indexofmaxpixel(hss, ch, CartesianIndices(size(hss.counts)[2:end]))

indexofmaxpixel(hss::HyperSpectrum) =
    indexofmaxpixel(hss, CartesianIndices(size(hss.counts)[2:end]))

indexofmaxpixel(hss::HyperSpectrum, cis::CartesianIndices) =
    map(i->indexofmaxpixel(hss,i,cis),Base.axes(hss,1))

"""
    Base.sum(hss::HyperSpectrum{T, N}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, mask::Union{BitArray, Array{Bool}}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, filt::Function) where {T<:Real, N}

Compute a sum spectrum for all or a subset of the pixels in `hss`.
"""

function Base.sum(hss::HyperSpectrum{T, N}, mask::Union{BitArray, Array{Bool}}) where {T<:Real, N}
    @assert size(mask)==size(hss) "$(size(mask)) ≠ $(size(hss))"
    RT = T isa Int ? Int64 : Float64
    vals = mapreduce(ci->hss.counts[:,ci], (a,b)->a .+= b, filter(ci->mask[ci], hss.index), init=zeros(RT, depth(hss)))
    props=copy(hss.properties)
    if haskey(hss.properties,:LiveTime)
        props[:LiveTime] = count(mask)*hss[:LiveTime]
    end
    return Spectrum(hss.energy, vals, props)
end

function Base.sum(hss::HyperSpectrum{T, N}) where {T<:Real, N}
    RT = T isa Int ? Int64 : Float64
    vals = mapreduce(ci->hss.counts[:,ci], (a,b)->a .+= b, hss.index, init=zeros(RT, depth(hss)))
    props=copy(hss.properties)
    if haskey(hss.properties,:LiveTime)
        props[:LiveTime] = length(hss)*hss[:LiveTime]
    end
    return Spectrum(hss.energy, vals, props)
end

function Base.sum(hss::HyperSpectrum{T, N}, filt::Function) where {T<:Real, N}
    RT = T isa Int ? Int64 : Float64
    vals = mapreduce(ci->hss.counts[:,ci], (a,b)->a .+= b, filter(ci->filt(hss,ci), hss.index), init=zeros(RT, depth(hss)))
    props=copy(hss.properties)
    if haskey(hss.properties,:LiveTime)
        props[:LiveTime] = count(ci->filt(hss,ci), hss.index) * hss[:LiveTime]
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
    ps = map(chs->plane(hss,chs,false),achs)
    maxval = maximum(map(p->maximum(p),ps))
    return map(p->Gray.(convert.(N0f8, p/maxval)),ps)
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
function roiimage(hss::HyperSpectrum, cxr::CharXRay, n=5)
    center = channel(energy(cxr), hss.energy)
    return roiimage(hss, center-n:center+n)
end

"""
    roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, n=5)

Create an array of Gray images representing the intensity in each of the CharXRay lines
in `cxrs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
function roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, n=5)
    achs = map( cxr->channel(energy(cxr), hss.energy)-n:channel(energy(cxr), hss.energy)+n, cxrs)
    return roiimages(hss, achs)
end
