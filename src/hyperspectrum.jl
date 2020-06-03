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

    :Elapse  # Total elapse time for map (so sig[:RealTime] ≈ sig[:Elapse]/length(sig)
    :Axes # Names for axes [ "Data", "Y", "X" ]
"""
struct Signal{T<:Real, N} <: AbstractArray{T, N}
    energy::EnergyScale
    properties::Dict{Symbol, Any}
    counts::Array{T, N}
    hash::UInt

    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, sz::Tuple{Int, Vararg{Int}}, ty::Type{<:Real}) =
         new{ty, length(sz)}(energy, properties, zeros(ty, sz), reduce(hash, hash.( (energy, properties, UInt(0x12347863)))))
    Signal(energy::EnergyScale, properties::Dict{Symbol,Any}, data::AbstractArray) =
        new{typeof(data[1]), ndims(data)}(energy, properties, data, mapreduce(hash,hash,(energy, properties, data)))
end

Base.show(io::IO, sig::Signal) =
    print(io,"Signal[$(get(sig.properties,:Name,"Unnamed")),$(sig.energy),$(size(sig)),$(size(sig.counts))]")

Base.getindex(sig::Signal, ci...) =
    getindex(sig.counts,ci...)
Base.getindex(sig::Signal, si::Symbol) =
    getindex(sig.properties, si)
Base.setindex!(sig::Signal, v::Real, ci::Int...) =
    setindex!(sig.counts, v, ci...)
Base.setindex!(sig::Signal, val::Any, si::Symbol) =
    setindex!(sig.properties, val, si)
Base.eltype(sig::Signal) = eltype(sig.counts)
Base.length(sig::Signal) = length(sig.counts)
Base.ndims(sig::Signal) = ndims(sig.counts)
Base.size(sig::Signal) = size(sig.counts)
Base.size(sig::Signal, n) = size(sig.counts, n)
Base.axes(sig::Signal) = Base.axes(sig.counts)
Base.axes(sig::Signal, n) = Base.axes(sig.counts, n)
Base.eachindex(sig::Signal) = eachindex(sig.counts)
Base.stride(sig::Signal, k) = stride(sig.counts, k)
Base.strides(sig::Signal) = strides(sig.counts)
depth(sig::Signal) = size(sig.counts,1)

NeXLCore.energy(sig::Signal, ch) = energy(ch, sig.energy)
channel(sig::Signal, energy) = channel(energy, sig.energy)

"""
    compressed(sig::Signal)

Returns a Signal with smaller or equal storage space to `sig` without losing any infomation.
"""
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

function Base.convert(ty::Type{<:AbstractFloat}, sig::Signal)
    res = convert(Array{ty,ndims(sig.counts)},sig.counts)
    return Signal(sig.energy, sig.properties, res)
end

"""
   plane(hss::Signal, chs::UnitRange, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the Signal.
"""
function plane(sig::Signal, chs::UnitRange, normalize=false)
    res = map(ci->sum(convert.(Float64, sig.counts[chs,ci])), CartesianIndices((size(sig)[2:end])))
    if normalize
        res/=maximum(res)
    end
    return res
end

"""
   plane(hss::Signal, ch::Int, normalize=false)

Sums a contiguous range of data planes into an Array. The dimension of the result is
one less than the dimension of the Signal.
"""
function plane(sig::Signal, ch::Int, normalize=false)
    res = sig.counts[ch,CartesianIndices((size(sig)[2:end]))]
    if normalize
        res/=maximum(res)
    end
    return res
end

"""
    maxpixel(sig::Signal)

Compute Bright's Max-Pixel derived signal.
"""
function maxpixel(sig::Signal)
    function maxi(i)
        tmp=zero(eltype(sig.counts))
        for idx in i:stride(sig.counts,2):length(sig.counts)
            tmp = max(tmp, sig.counts[idx])
        end
        return tmp
    end
    return map(i->maxi(i), Base.axes(sig.counts,1))
end

"""
    indexofmaxpixel(sig::Signal, ch::Int) # at channel `ch`
    indexofmaxpixel(sig::Signal) # all channels
    indexofmaxpixel(sig::Signal, ch::Int, cis::CartesianIndices)
    indexofmaxpixel(sig::Signal, cis::CartesianIndices)
    indexofmaxpixel(hs::HyperSpectrum, ch::Int) # at channel `ch`
    indexofmaxpixel(hs::HyperSpectrum) # all channels
    indexofmaxpixel(hs::HyperSpectrum, ch::Int, cis::CartesianIndices)
    indexofmaxpixel(hs::HyperSpectrum, cis::CartesianIndices)

Find the coordinates producing the maximum value in data[ch] or data[:] within 'cis' or full
spatial dimensions.
"""
function indexofmaxpixel(sig::Signal, ch::Int, cis::CartesianIndices)
    maxidx, max = cis[1], sig.counts[ch, cis[1]]
    for idx in cis
        if sig.counts[ch, idx] > max
            maxidx = idx
            max=sig.counts[ch, idx]
        end
    end
    return maxidx
end

indexofmaxpixel(sig::Signal, ch::Int) =
    indexofmaxpixel(sig, ch, CartesianIndices(size(sig.counts)[2:end]))

indexofmaxpixel(sig::Signal) =
    indexofmaxpixel(sig, CartesianIndices(size(sig.counts)[2:end]))

indexofmaxpixel(sig::Signal, cis::CartesianIndices) =
    map(i->indexofmaxpixel(sig,i,cis),Base.axes(sig,1))

"""
    sum(sig::Signal)

Computes the sum at each data index over all non-data axes.
"""
function Base.sum(sig::Signal)
    # Sum each plane (make sure result variable is large enough...)
    function sumi(i::Int)
        if eltype(sig.counts) in [ UInt8, UInt16, UInt32, Int8, Int16, Int32 ]
            tmp::Int64=0
            for idx in i:stride(sig.counts,2):length(sig.counts)
                @inbounds tmp += sig.counts[idx]
            end
            return tmp
        else
            tmpf::Float64=0.0
            for idx in i:stride(sig.counts,2):length(sig.counts)
                @inbounds tmpf += sig.counts[idx]
            end
            return tmpf
        end
    end
    return map(i->sumi(i), Base.axes(sig.counts,1))
end

"""
    sum(sig::Signal, filt::Function)

Produce a sum vector from those pixels for which `filt(hss, idx)==true` where `idx`
is an index over the dimensions 2:end.
"""
function Base.sum(sig::Signal, filt::Function)
    # Sum each plane (make sure result variable is large enough...)
    function sumi(i::Int, incld)
        if eltype(sig.counts) in [ UInt8, UInt16, UInt32, Int8, Int16, Int32 ]
            tmp::Int64=0
            for idx in incld
                tmp += sig.counts[i, idx]
            end
            return tmp
        else
            tmpf::Float64=0.0
            for idx in incld
                tmpf += sig.counts[i, idx]
            end
            return tmpf
        end
    end
    cis = CartesianIndices(( size(sig.counts)[2:end] ))
    include = filter(i->filt(sig, i), cis)
    return map(i->sumi(i, include), Base.axes(sig.counts,1))
end

"""
    roiimages(hss::Signal, achs::Vector{UnitRange{Int}})

Create an array of Gray images representing the intensity in each range of channels in
in `achs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
function roiimages(hss::Signal, achs::Vector{UnitRange{Int}})
    ps = map(chs->plane(hss,chs,false),achs)
    maxval = maximum(map(p->maximum(p),ps))
    return map(p->Gray.(convert.(N0f8, p/maxval)),ps)
end

"""
    HyperSpectrum

HyperSpectrum is a wrapper around Signal to facilitate access to the
the data as Spectrum objects.
"""
struct HyperSpectrum{T<:Real, N} <: AbstractArray{Spectrum{T}, N}
    signal::Signal{T}
    index::CartesianIndices # Spectrum indices
    properties::Dict{Symbol, Any}

    HyperSpectrum(sig::Signal) =
        new{typeof(sig.counts[1]),ndims(sig.counts)-1}(sig, CartesianIndices(size(sig.counts)[2:end]),sig.properties)
end

"""
    HyperSpectrum

Convert the Array{<:Real, N} perspective into a Array{Spectrum{<:Real}, N} perspective.

Special property tags:

    :Cartesian # The pixel index of a Spectrum extracted from a HyperSpectrum
"""
function ashyperspectrum(sig::Signal, name::AbstractString="Hyper-Spectrum")
    res=HyperSpectrum(sig)
    res[:Name]=name
    return res
end

Base.show(io::IO, hss::HyperSpectrum) =
    print(io,"HyperSpectrum[$(get(hss.signal.properties,:Name,"Unnamed")),$(hss.signal.energy),$(size(hss.index))]")

Base.eltype(hss::HyperSpectrum) = Spectrum
Base.length(hss::HyperSpectrum) = length(hss.index)
Base.ndims(hss::HyperSpectrum) = ndims(hss.index)
Base.size(hss::HyperSpectrum) = size(hss.index)
Base.size(hss::HyperSpectrum, n) = size(hss.index, n)
Base.axes(hss::HyperSpectrum) = Base.axes(hss.index)
Base.axes(hss::HyperSpectrum, n) = Base.axes(hss.index, n)
Base.eachindex(hss::HyperSpectrum) = hss.index
Base.stride(hss::HyperSpectrum, k) = stride(hss.index, k)
Base.strides(hss::HyperSpectrum) = strides(hss.index)
depth(hss::HyperSpectrum) = depth(hss.signal)
dose(hs::HyperSpectrum) = hs[:ProbeCurrent]*hs[:LiveTime]


function Base.getindex(hss::HyperSpectrum, idx...)::Spectrum
    so(i) = i*size(hss.signal.counts, 1)
    props = copy(hss.signal.properties)
    props[:Cartesian] = [ idx... ]
    return Spectrum(hss.signal.energy, hss.signal[:, idx...], props)
end

Base.getindex(hss::HyperSpectrum, i::Int)::Spectrum =
    getindex(hss, hss.index[i]...)
Base.getindex(hss::HyperSpectrum, sy::Symbol) =
    getindex(hss.signal.properties, sy)
Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol) =
    setindex!(hss.signal.properties, val, sy)

NeXLCore.energy(hs::HyperSpectrum, ch) = energy(ch, hs.signal.energy)
channel(hs::HyperSpectrum, energy) = channel(energy, hs.signal.energy)
rangeofenergies(hs::HyperSpectrum, ch) =
    ( energy(ch, hs.signal.energy), energy(ch+1, hs.signal.energy) )

"""
    sum(hss::HyperSpectrum)

Produce a sum spectrum.
"""
Base.sum(hss::HyperSpectrum)=
    Spectrum(hss.signal.energy, sum(hss.signal), copy(hss.signal.properties))

"""
    sum(sig::Signal, filt::Function)

Produce a sum Spectrum from those pixels for which filt(hss, idx)==true.
"""
Base.sum(hss::HyperSpectrum, filt::Function) =
    Spectrum(hss.signal.energy, sum(hss.signal, filt), copy(hss.signal.properties))

"""
    maxpixel(hss::HyperSpectrum)

Produce a maxpixel spectrum.
"""
maxpixel(hss::HyperSpectrum) =
    Spectrum(hss.signal.energy, maxpixel(hss.signal), copy(hss.signal.properties))

"""
    plane(hss::HyperSpectrum, chs::Union{Int,UnitRange{Int}}, norm=false) =

Extract as an Array the sum of the data in `chs`.
"""
plane(hss::HyperSpectrum, chs, norm=false) =
    plane(hss.signal, chs, norm)

"""
    roiimage(hss::Signal, chs::UnitRange)

Create a count map from the specified contiguous range of channels.
"""
roiimage(hss::HyperSpectrum, chs::UnitRange) =
    Gray.(convert.(N0f8, plane(hss, chs, true)))

"""
    roiimage(hss::Signal, cxr::CharXRay, n=5)

Create a count map for the specified characteristic X-ray.
"""
function roiimage(hss::HyperSpectrum, cxr::CharXRay, n=5)
    center = channel(energy(cxr), hss.signal.energy)
    return roiimage(hss, center-n:center+n)
end

"""
    roiimages(hss::HyperSpectrum, achs::Vector{UnitRange{Int}})

Create an array of Gray images representing the intensity in each range of channels in
in `achs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
roiimages(hss::HyperSpectrum, achs::Vector{UnitRange{Int}}) =
    roiimages(hss.signal, achs)

"""
    roiimages(hss::HyperSpectrum, cxrs::Vector{CharXRay}, n=5)

Create an array of Gray images representing the intensity in each of the CharXRay lines
in `cxrs`.  They are normalized such the the most intense pixel in any of them defines white.
"""
function roiimages(hss::HyperSpectrum, cxrs::Vector{CharXRay}, n=5)
    achs = map( cxr->channel(energy(cxr), hss.signal.energy)-n:channel(energy(cxr), hss.signal.energy)+n, cxrs)
    return roiimages(hss.signal, achs)
end

indexofmaxpixel(hs::HyperSpectrum, ch::Int, cis::CartesianIndices) =
    indexofmaxpixel(hs.signal, ch, cis)

indexofmaxpixel(hs::HyperSpectrum, ch::Int) =
    indexofmaxpixel(hs.signal, ch)

indexofmaxpixel(hs::HyperSpectrum, cis::CartesianIndices) =
    indexofmaxpixel(hs.signal, cis)

indexofmaxpixel(hs::HyperSpectrum) =
    indexofmaxpixel(hs.signal)
