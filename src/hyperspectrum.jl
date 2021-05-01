let hyperspectrumIndex = 0
    global hyperspectrumCounter() = (hyperspectrumIndex += 1)
end


"""
   HyperSpectrum(arr::Array{T<:Real}, energy::EnergyScale, props::Array{Symbol, Any})

   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}; axisnames = ( "X", "Y", "Z", "A", "B", "C" ), fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ), livetime=)
   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::AxisArray)
   HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, dims::NTuple{<:Integer}, depth::Int, type::Type{Real}; axisnames = ( "X", "Y", "Z", "A", "B", "C" ), fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )

The HyperSpectrum struct represents a multi-dimensional array of Spectrum objects.  The dimension of a HyperSpectrum
may be 1 for a traverse or a line-scan, 2 for a spectrum image or >2 for time-series of spectrum images or
multi-slice spectrum images.

The first constructor is used to create a HyperSpectrum from a raw Array of data.  The second to construct a HyperSpectrum
from another HyperSpectrum or an AxisArray.  The third from a description of the intended contents.

  * `axisnames`: A list of the names by which the axis can be referred
  * `fov`: The full width of the dimension in mm.

HyperSpectra differ from Array{Spectrum} in that the spectra in a HyperSpectrum must share properties like
:ProbeCurrent and :BeamEnergy.  However, each pixel can have it's own livetime.  HyperSpectrum objects can refer 
to line-scans (1D), spectrum images (2D), slice-and-view (3D), time sequenced images (3D), or higher dimension spectrum images.

Internally, HyperSpectrum reinterpretes an Array{T<:Real, N+1} as an Array{Spectrum{T<:Real},N-1}.

HyperSpectrum objects can be read from a RPL/RAW file (using `readrplraw(filenamebase::AbstractString)`) but
can be constructed from any Array{<:Real}.

HyperSpectrum objects can be indexed using the standard Julia array idioms including a single integer index or a
CartesianIndex.  For example, to iterate over every spectrum in a HyperSpectrum

    % Construct a 21 row × 19 column spectrum image with 2048 channels of [0,255].
    hs = HyperSpectrum{es, props, (21,19), 2048, UInt8}
    for idx in eachindex(hs)
        spec = hs[idx]   % get a Spectrum representing the 2048 channels of data at idx
        spec[22] = 1     % Set the 22nd channel to 1
    end
    % Indices into a HyperSpectrum are (row, column)
"""
struct HyperSpectrum{T<:Real,N} <: AbstractArray{Spectrum{T},N}
    counts::AxisArray{T} # Organized as (ch, r, c) or (ch, y, x)
    index::CartesianIndices # Spectrum indices
    energy::EnergyScale
    livetime::AbstractArray{Float64} # A matrix containing the pixel-by-pixel livetime
    properties::Dict{Symbol,Any}
    stagemap::Type{<:StageMapping} # maps pixel coordinates into stage coordinates

    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}; 
        axisnames = ( "Y", "X", "Z", "A", "B", "C" ), #
        fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ), #
        stagemap::Type{<:StageMapping}=DefaultStageMapping,
        livetime=get(props, :LiveTime, 1.0) * ones(Float64, size(arr)[2:end]...)
    )
        axes = [ Axis{:Channel}(1:size(arr,1)), #
                ( Axis{Symbol(axisnames[i-1])}(fov[i-1]*(-0.5:1.0/(size(arr,i)-1):0.5)) for i in 2:ndims(arr) )... ]
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        new{eltype(arr),ndims(arr) - 1}(
            AxisArray(arr, axes...),
            CartesianIndices(size(arr)[2:end]),
            energy,
            livetime,
            props,
            stagemap,
        )
    end
    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::AxisArray; stagemap::Type{<:StageMapping}=DefaultStageMapping, livetime=get(props, :LiveTime, 1.0)  * ones(Float64, size(arr[2:end])...))
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        new{eltype(arr),ndims(arr) - 1}(
            arr,
            CartesianIndices(size(arr)[2:end]),
            energy,
            livetime,
            props,
            stagemap
        )
    end
    function HyperSpectrum(
        energy::EnergyScale,
        props::Dict{Symbol,Any},
        dims::Tuple,
        depth::Int,
        type::Type{<:Real};
        axisnames = ( "X", "Y", "Z", "A", "B", "C" ), 
        fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ),
        stagemap::Type{<:StageMapping}=DefaultStageMapping
    )
        axes = [ Axis{:Channel}(1:depth), #
                ( Axis{Symbol(axisnames[i])}(fov[i]*(-0.5:1.0/(dims[i]-1):0.5)) for i in eachindex(dims) )... ]
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        new{type, length(dims)}(
            AxisArray(zeros(type, depth, dims...), axes...),
            CartesianIndices(dims),
            energy,
            get(props,:LiveTime,1.0) * ones(Float64, dims...),
            props,
            stagemap
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$(eltype(hss.counts)),$(axisnames(hss))}[$(name(hss)), $(hss.energy), $(depth(hss)) ch]",
    )
end

function Base.show(io::IO, ::MIME"text/html", hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "<p>$sz HyperSpectrum{$(eltype(hss.counts)),$(axisnames(hss))}[$(name(hss))), $(hss.energy), $(depth(hss)) ch]</p>",
    )
end

function Base.show(io::IO, hss::HyperSpectrum)
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$(eltype(hss.counts)),$(axisnames(hss))}[$(name(hss)), $(hss.energy), $(depth(hss)) ch]",
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
Base.stride(hss::HyperSpectrum, k::Int) = stride(hss.counts.data, k + 1)
Base.strides(hss::HyperSpectrum) = strides(hss.counts.data)[2:end]
depth(hss::HyperSpectrum) = size(hss.counts, 1)
NeXLCore.name(hss::HyperSpectrum) = get(hss.properties, :Name, "HyperSpectrum")
AxisArrays.axisnames(hss::HyperSpectrum) = AxisArrays.axisnames(hss.counts)[2:end]
AxisArrays.axisvalues(hss::HyperSpectrum) = AxisArrays.axisvalues(hss.counts)[2:end]
axisname(hss::HyperSpectrum, ax::Int) = AxisArrays.axisnames(hss.counts)[ax+1]
axisvalue(hss::HyperSpectrum, ax::Int, j::Int) = AxisArrays.axisvalues(hss.counts)[ax+1][j]
axisrange(hss::HyperSpectrum, ax::Int) = AxisArrays.axisvalues(hss.counts)[ax+1]

"""
    coordinate(hss::HyperSpectrum, idx::Tuple{<:Int})

Computes the stage coordinate centering the pixel specified by `idx` using the `StageMapping` `hss.stagemap`.
"""
function coordinate(hss::HyperSpectrum, idx::Tuple)
    av, an = axisvalues(hss.counts), axisnames(hss.counts)
    img_coord = Dict{Symbol,Float64}(  an[i+1] => av[i+1][ii] for (i,ii) in enumerate(idx) )
    stage_coord = get(hss.properties, :StagePosition, Dict{Symbol,Float64}(:X=>0.0, :Y=>0.0))
    th = get(hss.properties, :ImageRotation, 0.0)
    return image2stage(hss.stagemap, stage_coord, img_coord, th)
end
coordinate(hss::HyperSpectrum, ci::CartesianIndex) = #
    coordinate(hss, ci.I)

"""
    dose(hss::HyperSpectrum)
    dose(hss::HyperSpectrum, idx::CartesianIndex)
    dose(hss::HyperSpectrum, idx...)

Returns the product of the probe current and the live-time on a per pixel basis.
"""
dose(hss::HyperSpectrum) = hss[:ProbeCurrent] * mean(hss.livetime)
dose(hss::HyperSpectrum, idx::CartesianIndex) = hss[:ProbeCurrent] * hss.livetime[idx]
dose(hss::HyperSpectrum, idx::Int...)  = hss[:ProbeCurrent] * hss.livetime[idx...]


"""
    livetime!(hss::HyperSpectrum,idx::CartesianIndex, lt::AbstractFloat)
    livetime!(hss::HyperSpectrum{T,N}, idx::NTuple{Int,N}, lt::AbstractFloat)
    livetime!(hss::HyperSpectrum{T,N}, lt::AbstractFloat) # All pixels to lt

Set the livetime on a per pixel basis.
"""
function livetime!(hss::HyperSpectrum, idx::CartesianIndex, lt::AbstractFloat) 
    hss.livetime[idx] = lt
end
function livetime!(hss::HyperSpectrum{T,N}, idx::NTuple{N, Int}, lt::AbstractFloat)  where {T<:Real,N}  
    hss.livetime[idx...] = lt
end
function livetime!(hss::HyperSpectrum, lt::AbstractFloat)  
    hss.livetime .= lt
end

livetime(hss::HyperSpectrum, idx::CartesianIndex) = hss.livetime[idx]
livetime(hss::HyperSpectrum{T,N}, idx::NTuple{N, Int}) where {T<:Real,N} = hss.livetime[idx...]
livetime(hss::HyperSpectrum, idx::Int...) = hss.livetime[idx...]


function Base.getindex(
    hss::HyperSpectrum{T,N},
    ci::CartesianIndex
)::Spectrum{T} where { T<: Real, N} 
    props = copy(hss.properties)
    props[:Cartesian] = ci
    props[:Name] = "$(name(hss))[$(join(["$i" for i in ci.I],","))]"
    props[:StageCoordinate]=coordinate(hss, ci)
    props[:LiveTime] = hss.livetime[ci]
    return Spectrum(hss.energy, hss.counts[:, ci], props)
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    idx::Vararg{Int,N},
)::Spectrum{T} where {T<:Real,N}
    getindex(hss, CartesianIndex(idx...))
end
function Base.getindex(hss::HyperSpectrum{T,N}, idx::Int)::Spectrum{T} where {T<:Real,N}
    getindex(hss, hss.index[idx])
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    vidx::Vector
)::Vector{Spectrum{T}} where { T<: Real, N} 
    return [ hss[i] for i in vidx ]
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    rng::AbstractRange{<:Integer},
)::Vector{Spectrum{T}} where { T<: Real, N} 
    return [ hss[i] for i in rng ]
end
function Base.getindex(
    hss::HyperSpectrum{T,N},
    idx...
)::HyperSpectrum{T,N} where {T<:Real, N}
    props = copy(hss.properties)
    props[:Name] = "$(name(hss))[$idx]"
    return HyperSpectrum(hss.energy, props, hss.counts[:, idx...], livetime=hss.livetime[idx...])
end

function Base.getindex(
    hss::HyperSpectrum,
    cxr::CharXRay
) 
    return roiimage(hss, cxr)
end
function Base.getindex(
    hss::HyperSpectrum,
    cxrs::Vector{CharXRay}
) 
    return roiimages(hss, cxrs)
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
    HyperSpectrum(hss.energy, hss.properties, hss.counts[:, ranges...], livetime=hss.livetime[ranges...])
end

function Base.getindex(hss::HyperSpectrum, sy::Symbol)
    @assert sy!=:LiveTime "Use the livetime(...) function to get the :LiveTime property for hyperspectra."
    getindex(hss.properties, sy)
end
function Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol)
    @assert sy!=:LiveTime "Use the livetime!(...) function to set the :LiveTime property for hyperspectra."
    setindex!(hss.properties, val, sy)
end

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
counts(hss::HyperSpectrum, ci::CartesianIndex, ch::Int) = hss.counts[ch, ci]
counts(hss::HyperSpectrum, ci::CartesianIndex) = view(hss.counts[:, ci])

"""
    compress(hss::HyperSpectrum)

Returns a HyperSpectrum with smaller or equal storage space to `hss` without losing or truncating any counts 
(note: AbstractFloat compresses to Float32 with loss of precision).  Can change the storage type and/or 
reduce the depth of hss.
"""
function compress(hss::HyperSpectrum{T,N}) where {T<:Real,N}
    data = hss.counts
    (minval, maxval) = extrema(data)
    maxi = mapreduce(max, CartesianIndices(hss)) do ci
        something(findlast(c -> c != 0.0, data[:, ci]), 0)
    end
    data = maxi < depth(hss) ? data[1:maxi, axes(hss.counts)[2:end]...] : data
    if eltype(T) in (Int16, Int32, Int64)
        for newtype in filter(ty -> sizeof(T) > sizeof(ty), (Int8, Int16, Int32))
            if minval >= typemin(newtype) && maxval <= typemax(newtype)
                return HyperSpectrum(hss.energy, copy(hss.properties), newtype.(data), livetime=hss.livetime)
            end
        end
    elseif T in (UInt16, UInt32, UInt64)
        for newtype in filter(ty -> sizeof(T) > sizeof(ty), (UInt8, UInt16, UInt32))
            if maxval <= typemax(newtype)
                return HyperSpectrum(hss.energy, copy(hss.properties), newtype.(data), livetime=hss.livetime)
            end
        end
    elseif T isa Float64
        if sizeof(T) < sizeof(Float32) &&
           minval >= typemin(Float32) &&
           maxval <= typemax(Float32)
            return HyperSpectrum(hss.energy, copy(hss.properties), Float32.(data), livetime=hss.livetime)
        end
    end
    return maxi < depth(hss) ? HyperSpectrum(hss.energy, copy(hss.properties), data, livetime=hss.livetime) : hss
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
    res = map(ci -> sum(Float64.(hss.counts[chs, ci])), hss.index)
    return normalize ? res / maximum(res) : res
end

"""
   plane(hss::HyperSpectrum, ch::Int, normalize=false)

Extracts a single plane from a HyperSpectrum. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, ch::Int, normalize = false)
    res = hss.counts[ch, hss.index]
    return normalize ? res / maximum(res) : res
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
maxpixel(hss::HyperSpectrum, mask::BitArray) =
    maxpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function maxpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
)::Spectrum{T} where {T<:Real,N}
    data, res = hss.counts, zeros(T, depth(hss))
    for ci in filter(ci -> filt(ci), cis)
        res .= max.(res, data[:, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res, props)
end
function maxpixel(
    hss::HyperSpectrum{T,N},
)::Spectrum{T} where {T<:Real,N}
    res = map(Base.OneTo(size(hss.counts, 1))) do ch
        maximum(view(hss.counts, ch:size(hss.counts, 1):prod(size(hss.counts))))
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res, props)
end

function minpixel(
    hss::HyperSpectrum{T,N},
)::Spectrum{T} where {T<:Real,N}
    res = map(Base.OneTo(size(hss.counts, 1))) do ch
        minimum(view(hss.counts, ch:size(hss.counts, 1):prod(size(hss.counts))))
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res, props)
end
minpixel(hss::HyperSpectrum, mask::BitArray) =
    minpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function minpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real,N}
    data, res = hss.counts, fill(typemax(T), (depth(hss),))
    for ci in filter(ci -> filt(ci), cis), ch in Base.OneTo(depth(hss))
        res[ch] = min(res[ch], data[ch, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res, props)
end
avgpixel(hss::HyperSpectrum, mask::BitArray) =
    avgpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function avgpixel(
    hss::HyperSpectrum{T,N},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real,N}
    data, res, cx = hss.counts, zeros(Float64, depth(hss)), 0
    for ci in filter(ci -> filt(ci), cis)
        for ch in Base.OneTo(depth(hss))
            res[ch] += data[ch, ci]
        end
        cx += 1
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res ./ max(1.0, cx), props)
end
function avgpixel(
    hss::HyperSpectrum{T,N}
) where {T<:Real,N}
    res = map(Base.OneTo(size(hss.counts, 1))) do ch
        mean(view(hss.counts, ch:size(hss.counts, 1):prod(size(hss.counts))))
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    return Spectrum(hss.energy, res, props)
end

function sumcounts(hss::HyperSpectrum, cis::CartesianIndices = CartesianIndices(hss))
    return [sum(Int.(hss.counts[:, ci])) for ci in cis]
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
    data = hss.counts
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
    data = hss.counts
    for ci in cis, i in filter(i -> data[i, ci] > res[i], eachindex(res))
        res[i] = data[i, ci]
        cix[i] = ci
    end
    return cix
end

"""
    Base.sum(hss::HyperSpectrum{T, N}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, mask::Union{BitArray, Array{Bool}}) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, func::Function) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, cis::CartesianIndices) where {T<:Real, N}

Compute a sum spectrum for all or a subset of the pixels in `hss`.

Where func(hss::HyperSpectrum, ci::CartesianIndex)::Bool
"""
function Base.sum(
    hss::HyperSpectrum{T,N},
    mask::Union{BitArray,Array{Bool}};
    name = nothing
) where {T<:Real,N}
    @assert size(mask) == size(hss) "Mask size[$(size(mask))] ≠ Hyperspectrum size[$(size(hss))]"
    data = hss.counts
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in hss.index
        if mask[ci]
            res .+= data[:, ci]
            lt += hss.livetime[ci]
        end
    end
    props = copy(hss.properties)
    props[:Name] = something(name, props[:Name])
    return Spectrum(hss.energy, res, props)
end
function Base.sum(
    hss::HyperSpectrum{T,N},
    filt::Function;
    name = nothing
) where {T<:Real,N}
    data = hss.counts
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in hss.index
        if filt(hss,ci)
            res .+= data[:, ci]
            lt += hss.livetime[ci]
        end
    end
    props = copy(hss.properties)
    props[:Name] = something(name, props[:Name])
    return Spectrum(hss.energy, res, props)
end
function Base.sum(hss::HyperSpectrum{T,N}) where {T<:Real,N}
    data = hss.counts
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in hss.index
        res .+= data[:, ci]
        lt += hss.livetime[ci]
    end
    props = copy(hss.properties)
    props[:Name] = "Sum[$hss]"
    return Spectrum(hss.energy, res, props)
end
function Base.sum(hss::HyperSpectrum{T,N}, cis::CartesianIndices) where {T<:Real,N}
    data = hss.counts
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in cis
        res .+= data[:, ci]
        lt += hss.livetime[ci]
    end
    props = copy(hss.properties)
    props[:Name] = "Sum[$hss, $cis]"
    return Spectrum(hss.energy, res, props)
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
