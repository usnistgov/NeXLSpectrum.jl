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
    counts::Array{T,N}

    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, sz::Tuple{Int, Vararg{Int}}, ty::Type{<:Real}) =
         new{ty, Int}(energy, properties, zeros(ty, sz))
    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, data::AbstractArray) =
        new{typeof(data[1]), Int}(energy, properties, data)

end

Base.show(io::IOStream, sig::Signal) =
    print(io,"Signal[$(sig.energy),$(size(sig.counts))]")

Base.getindex(sig::Signal, ci::CartesianIndex) =
    sig.counts[ci]
Base.setindex!(sig::Signal, v::Real, ci::CartesianIndex) =
    setindex(sig.counts, v, ci)
Base.setindex!(sig::Signal, v::Real, ci::Tuple{Integer, Vararg{Integer}}) =
    setindex(sig.counts, v, CartesianIndex(ci))
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

"""
    writecounts(ios::IOStream, hs::Signal)

Write the counts data directly in binary to the specified stream.
"""
writecounts(ios::IOStream, hs::Signal) =
    write(ios, hs.counts)

"""
    readcounts(ios::IOStream, hs::Signal)

Read counts data directly in binary from the specified stream.
"""
readcounts(ios::IOStream, T::Type{<:Real}, size::Tuple{Integer, Vararg{Integer}}) =
    read!(ios, Array{T}(undef, size))


"""
    readraw(ios::IOStream, T::Type{<:Real}, size::NTuple{Int, N}, energy::EnergyScale, props::Dict{Symbol,Any}) =

Construct a full Signal from the counts data in binary from the specified stream.
"""
readraw(ios::IOStream, T::Type{<:Real}, size::Tuple{Integer, Vararg{Integer}}, energy::EnergyScale, props::Dict{Symbol,Any}) =
    Signal(energy, props, readcounts(ios, T, size))


"""
    HyperSpectrum

HyperSpectrum is a wrapper around Signal to facilitate access of the
the data as individual Spectrum objects.
"""
struct HyperSpectrum{T<:Real, N} <: AbstractArray{Spectrum{T}, N}
    hs::Signal{T, N}
    index::CartesianIndices # Spectrum indices

    HyperSpectrum(hs::Signal) =
        new{typeof(hs.counts[1]),Int}(hs, CartesianIndices(size(hs.counts)[2:end]))
end

"""
    HyperSpectrum

Convert the Array{<:Real, N} perspective into a Array{Spectrum{<:Real}, N} perspective.

Special property tags:

    :Cartesian # The pixel index of a Spectrum extracted from a HyperSpectrum
"""
ashyperspectrum(hs::Signal) = HyperSpectrum(hs)


Base.eltype(hss::HyperSpectrum) = Spectrum
Base.length(hss::HyperSpectrum) = length(hs.index)
Base.ndims(hss::HyperSpectrum) = ndims(hs.index)
Base.size(hss::HyperSpectrum) = size(hs.index)
Base.size(hss::HyperSpectrum, n) = size(hs.index, n)
Base.axes(hss::HyperSpectrum) = axes(hs.index)
Base.axes(hss::HyperSpectrum, n) = axes(hs.index, n)
Base.eachindex(hss::HyperSpectrum) = hs.index
Base.stride(hss::HyperSpectrum, k) = stride(hs.index, k)
Base.strides(hss::HyperSpectrum) = strides(hs.index)

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

"""
   plane(hss::HyperSpectrum, chs::UnitRange)

Sum a consecutive series of channel planes into a single Array representing an image plane.
"""
function plane(hs::Signal, chs::UnitRange)
    res = Array{Float64}(undef, size(hs)[2:end])
    for idx in eachindex(res)
        res[i] = sum(hs.counts[chs, idx])
    end
    return res
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
