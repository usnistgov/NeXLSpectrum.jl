using Images

"""
    Signal

The multidimensional equivalent of Spectrum.  A megapixel spectrum image
might be constructed as Signal(energy, props, (4096, 1024, 2048)) where
there are 4096 channels, 1024 rows and 2048 columns.  HyperSpectra may be 1, 2,
.. N dimensional but since they reside in memory, there are practical limits.

A type for data sets containing multiple closely related spectra.  A Signal is
slighly more restricted than an Array{Spectrum,N} because all the Spectra in a Signal are
assumed to have certain properties in common -  the EnergyScale and a set of common Spectrum
properties.  You can also specify an efficient packing of the data by using
UInt8, UInt16, ... etc as required to hold the data.

Signal maintains the data as an abstract array but provides functions to
extract individual points or sets of points as Spectrum.

Special property tags:

    :Elapse  # Total elapse time for map (so hs[:RealTime] ≈ hs[:Elapse]/length(hs)
"""
struct Signal{T<:Real, N} <: AbstractArray{T, N}
    energy::EnergyScale
    properties::Dict{Symbol, Any}
    counts::Array{T, N}
    hash::UInt

    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, sz::Tuple{Int, Vararg{Int}}, ty::Type{<:Real}) =
         new{ty, Int}(energy, properties, zeros(ty, sz), reduce(hash, hash.(energy, properties, UInt(0x12347863))))
    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, data::AbstractArray) =
        new{typeof(data[1]), ndims(data)}(energy, properties, data, mapreduce(hash,hash,(energy, properties, data)))
end

Base.show(io::IO, sig::Signal) =
    print(io,"Signal[$(sig.energy),$(size(sig.counts))]")

Base.getindex(sig::Signal, ci::CartesianIndex) =
    sig.counts[ci]
Base.getindex(sig::Signal, ci...) =
    sig.counts[ci...]
Base.setindex!(sig::Signal, v::Real, ci::CartesianIndex) =
    setindex(sig.counts, v, ci)
Base.setindex!(sig::Signal, v::Real, ci...) =
    setindex(sig.counts, v, ci...)
Base.eltype(sig::Signal) = eltype(sig.counts)
Base.length(sig::Signal) = length(sig.counts)
Base.ndims(sig::Signal) = ndims(sig.counts)
Base.size(sig::Signal) = size(sig.counts)
Base.size(sig::Signal, n) = size(sig.counts, n)
Base.axes(sig::Signal) = axes(sig.counts)
Base.axes(sig::Signal, n) = axes(sig.counts, n)
Base.eachindex(sig::Signal) = eachindex(sig.counts)
Base.stride(sig::Signal, k) = stride(sig.counts, k)
Base.strides(sig::Signal) = strides(sig.counts)

NeXLCore.energy(sig::Signal, ch) = energy(ch, sig.energy)
channel(sig::Signal, energy) = channel(energy, sig.energy)

function asimage(sig::Signal, ch)
    if typeof(hs.counts[1])==UInt8
        return normedview(N0f8,raw.counts[ch,:,:])
    elseif typeof(hs.counts[1])==UInt16
        return normedview(N0f16,raw.counts[ch,:,:])
    end
end

function compressed(sig::Signal)
    maxval = maximum(sig.counts)
    if (typeof(sig.counts[1]) in ( UInt16, UInt32 )) && (maxval <= 255)
        res = convert(Array{UInt8,ndims(sig.counts)},sig.counts)
        return Signal(sig.energy, sig.properties, res)
    elseif (typeof(sig.counts[1]) isa UInt32 ) && (maxval <= 65535)
        res = convert(Array{UInt16,ndims(sig.counts)},sig.counts)
        return Signal(sig.energy, sig.properties, res)
    elseif (typeof(sig.counts[1]) in ( Int16, Int32 )) && (maxval <= 127)
        res = convert(Array{Int8,ndims(sig.counts)},sig.counts)
        return Signal(sig.energy, sig.properties, res)
    elseif (typeof(sig.counts[1]) isa Int32 ) && (maxval <= 32767)
        res = convert(Array{Int16,ndims(sig.counts)},sig.counts)
        return Signal(sig.energy, sig.properties, res)
    elseif typeof(sig.counts[1]) isa Float64
        res = convert(Array{Float32,ndims(sig.counts)},sig.counts)
        return Signal(sig.energy, sig.properties, res)
    end
    return sig
end

"""
   plane(hss::Signal, chs::UnitRange, normalize=false)

Sum a consecutive series of channel planes into a single Array representing an image plane.
"""
function plane(sig::Signal, chs::UnitRange, normalize=false)
    res = Array{Float64}(undef, size(hs)[2:end])
    for idx in eachindex(res)
        res[i] = sum(hs.counts[chs, idx])
    end
    if normalize
        res/=maximum(res)
    end
    return res
end

"""
    HyperSpectrum

HyperSpectrum is a wrapper around Signal to facilitate access to the
the data as Spectrum objects.
"""
struct HyperSpectrum{T<:Real, N} <: AbstractArray{Spectrum{T}, N}
    hs::Signal{T}
    index::CartesianIndices # Spectrum indices

    HyperSpectrum(sig::Signal) =
        new{typeof(sig.counts[1]),ndims(sig.counts)-1}(sig, CartesianIndices(size(sig.counts)[2:end]))
end

"""
    HyperSpectrum

Convert the Array{<:Real, N} perspective into a Array{Spectrum{<:Real}, N} perspective.

Special property tags:

    :Cartesian # The pixel index of a Spectrum extracted from a HyperSpectrum
"""
ashyperspectrum(sig::Signal) = HyperSpectrum(sig)


Base.eltype(hss::HyperSpectrum) = Spectrum
Base.length(hss::HyperSpectrum) = length(hss.index)
Base.ndims(hss::HyperSpectrum) = ndims(hss.index)
Base.size(hss::HyperSpectrum) = size(hss.index)
Base.size(hss::HyperSpectrum, n) = size(hss.index, n)
Base.axes(hss::HyperSpectrum) = axes(hss.index)
Base.axes(hss::HyperSpectrum, n) = axes(hss.index, n)
Base.eachindex(hss::HyperSpectrum) = hss.index
Base.stride(hss::HyperSpectrum, k) = stride(hss.index, k)
Base.strides(hss::HyperSpectrum) = strides(hss.index)

function Base.getindex(hss::HyperSpectrum, ci::CartesianIndex)::Spectrum
    so(i) = i*size(hss.hs.counts, 1)
    props = copy(hss.hs.props)
    props[:Cartesian] = [ ci.I ]
    return Spectrum(hss.hs.energy, props, hss.hs.counts[:, ci])
end

Base.getindex(hss::HyperSpectrum, i::Int)::Spectrum =
    getindex(hss::HyperSpectrum, hss.index[i])

"""
    sum(hss::HyperSpectrum, filt::Function)

Produce a sum spectrum from those pixels at idx for which filt(hss, idx)==true.
"""
function Base.sum(hss::HyperSpectrum, filt::Function)
    include = map(i->filt(hss, i), hss.index)
    a = copy(hss.counts[:, include[1]])
    cis = [ include[1].I ]
    for incl in include[2:end]
        a .+= hss.counts[:, incl]
        push!(cis, incl.I)
    end
    props = copy(hss.properties)
    if haskey(props, :LiveTime)
        props[:LiveTime] = length(include)*props[:LiveTime]
    end
    return Spectrum(hss.energy, props, a)
end


plane(hs::HyperSpectrum,chs::UnitRange) =
    plane(hss.sig,chs)
"""
    countmap(hss::Signal, chs::UnitRange)

Create a count map from the specified contiguous range of channels.
"""
function countmap(hss::HyperSpectrum, chs::UnitRange)
    pl = plane(hss, chs)
    return Gray.(pl ./ max(pl))
end
