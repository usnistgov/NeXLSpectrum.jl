using Images

"""
    HyperSpectrum

The multidimensional equivalent of Spectrum.  A megapixel spectrum image
might be constructed as HyperSpectrum(energy, props, (4096, 1024, 2048)) where
there are 4096 channels, 1024 rows and 2048 columns.  HyperSpectra may be 1, 2,
.. N dimensional but since they reside in memory, there are practical limits.

A type for data sets containing multiple closely related spectra.  A HyperSpectrum is
slighly more restricted than an Array{Spectrum,N} because all the Spectra in a HyperSpectrum are
assumed to have certain properties in common -  the EnergyScale and a set of common Spectrum
properties.  You can also specify an efficient packing of the data by using
UInt8, UInt16, ... etc as required to hold the data.

HyperSpectrum maintains the data as an abstract array but provides functions to
extract individual points or sets of points as Spectrum.

Special property tags:

    :Elapse  # Total elapse time for map (so hs[:RealTime] ≈ hs[:Elapse]/length(hs)
"""
struct HyperSpectrum{T<:Real} <: AbstractArray{T, N}
    energy::EnergyScale
    properties::Dict{Symbol, Any}
    counts::Array{T, N}

    HyperSpectrum(energy::EnergeScale, properties::Dict{Symbol,Any}, sz::NTuple{N,Int}, T::Type{<:Real}=Int16) =
        new(energy, properties, Array{T}(zero(T),sz))
end

getindex(hs::HyperSpectrum, ci::CartesianIndex) =
    hs.counts[ci]
getindex(hs::HyperSpectrum, ci::Tuple{Vararg{Integer, N}}) =
    getindex(hsd.counts, CartesianIndex(ci))
setindex!(hs::HyperSpectrum, v::Real, ci::CartesianIndex) =
    setindex(hs.counts, v, ci)
setindex!(hs::HyperSpectrum, v::Real, ci::Tuple{Vararg{Integer, N}}) =
    setindex(hs.counts, v, CartesianIndex(ci))
eltype(hs::HyperSpectrum) = eltype(hs.counts)
length(hs::HyperSpectrum) = length(hs.counts)
ndims(hs::HyperSpectrum) = ndims(hs.counts)
size(hs::HyperSpectrum) = size(hs.counts)
size(hs::HyperSpectrum, n) = size(hs.counts, n)
axes(hs::HyperSpectrum) = axes(hs.counts)
axes(hs::HyperSpectrum, n) = axes(hs.counts, n)
eachindex(hs::HyperSpectrum) = eachindex(hs.counts)
stride(hs::HyperSpectrum, k) = stride(hs.counts, k)
stides(hs::HyperSpectrum) = strides(hs.counts)

energy(hs::HyperSpectrum, ch) = energy(ch, hs.energy)
channel(hs::HyperSpectrum, energy) = channel(energy, hs.energy)

"""
    HyperSpectrumS

Convert the Array{<:Real, N} perspective into a Array{Spectrum{<:Real}, N} perspective.

Special property tags:

    :Elapse  # Total elapse time for map (so hs[:RealTime] ≈ hs[:Elapse]/length(hs)
    :Cartesian # The pixel index of a Spectrum extracted from a HyperSpectrumS
"""
asspectra(hs::HyperSpectrum) = HyperSpectrumS(hs)

"""
    writecounts(ios::IOStream, hs::HyperSpectrum)

Write the counts data directly in binary to the specified stream.
"""
writecounts(ios::IOStream, hs::HyperSpectrum) =
    write(ios, hs.counts)

"""
    readcounts(ios::IOStream, hs::HyperSpectrum)

Read counts data directly in binary from the specified stream.
"""
readcounts(ios::IOStream, T::Type{<:Real}, size::NTuple{Int, N}) =
    read(ios, Array{T}(undef, size))


"""
    readraw(ios::IOStream, T::Type{<:Real}, size::NTuple{Int, N}, energy::EnergyScale, props::Dict{Symbol,Any}) =

Construct a full HyperSpectrum from the counts data in binary from the specified stream.
"""
readraw(ios::IOStream, T::Type{<:Real}, size::NTuple{Int, N}, energy::EnergyScale, props::Dict{Symbol,Any}) =
    HyperSpectrum(energy, props, readcounts(ios, T, size))


"""
    HyperSpectrumS

HyperSpectrumS is a wrapper around HyperSpectrum to facilitate access of the
the data as individual Spectrum objects.
"""

struct HyperSpectrumS <: AbstractArray{Spectrum, N}
    hs::HyperSpectrum
    index::CartesianIndices # Spectrum indices

    HyperSpectrumS(hs::HyperSpectrum) = new(hs,CartesianIndices(size(hs)[2:end]))
end

eltype(hss::HyperSpectrumS) = Spectrum
length(hss::HyperSpectrumS) = length(hs.index)
ndims(hss::HyperSpectrumS) = ndims(hs.index)
size(hss::HyperSpectrumS) = size(hs.index)
size(hss::HyperSpectrumS, n) = size(hs.index, n)
axes(hss::HyperSpectrumS) = axes(hs.index)
axes(hss::HyperSpectrumS, n) = axes(hs.index, n)
eachindex(hss::HyperSpectrumS) = hs.index
stride(hss::HyperSpectrumS, k) = stride(hs.index, k)
stides(hss::HyperSpectrumS) = strides(hs.index)

function getindex(hss::HyperSpectrumS, ci::CartesianIndex)::Spectrum
    so(i) = i*size(hss.hs.counts, 1)
    props = copy(hss.hs.props)
    props[:Cartesian] = [ ci.I ]
    return Spectrum(hss.hs.energy, props, hss.hs.counts[:, ci])
end

getindex(hss::HyperSpectrumS, i::Int)::Spectrum =
    getindex(hss::HyperSpectrumS, hss.index[i])

"""
    sum(hss::HyperSpectrumS, filt::Function)

Produce a sum spectrum from those pixels at idx for which filt(hss, idx)==true.
"""
function sum(hss::HyperSpectrumS, filt::Function)
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
   plane(hss::HyperSpectrumS, chs::UnitRange)

Sum a consecutive series of channel planes into a single Array representing an image plane.
"""
function plane(hss::HyperSpectrumS, chs::UnitRange)
    res = Array{Float64}(undef, size(hss.index))
    for idx in hss.index
        res[i] = sum(hss.hs.counts[chs, idx])
    end
    return res
end

"""
    countmap(hss::HyperSpectrum, chs::UnitRange)

Create a count map from the specified contiguous range of channels.
"""
function countmap(hss::HyperSpectrum, chs::UnitRange)
    pl = plane(hss, chs)
    return Gray.(pl ./ max(pl))
end
