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
struct HyperSpectrum{T<:Real, N, NP} <: AbstractArray{Spectrum{T}, N}
    counts::AxisArray{T, NP, Array{T,NP}} # Organized as (ch, r, c) or (ch, y, x)
    energy::EnergyScale
    livetime::Array{Float64, N} # A matrix containing the pixel-by-pixel livetime
    properties::Dict{Symbol,Any}
    stagemap::Type{<:StageMapping} # maps pixel coordinates into stage coordinates
    hash::UInt

    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::Array{<:Real}; 
        axisnames = ( "Y", "X", "Z", "A", "B", "C" ), #
        fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ), #
        offset = ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ), #
        stagemap::Type{<:StageMapping}=DefaultStageMapping, #
        livetime=fill(get(props, :LiveTime, 1.0), size(arr)[2:end]...)
    )
        axes = [ Axis{:Channel}(1:size(arr,1)), #
                ( Axis{Symbol(axisnames[i-1])}(offset[i-1]:fov[i-1]/(size(arr,i)-1):offset[i-1]+fov[i-1]) for i in 2:ndims(arr) )... ]
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        counts = AxisArray(arr, axes...)
        # @show AxisArrays.axes(counts)
        new{eltype(arr), ndims(arr) - 1, ndims(arr)}(
            counts,
            energy,
            livetime,
            props,
            stagemap,
            hash(counts, hash(energy, hash(livetime, hash(props, hash(stagemap)))))
        )
    end
    function HyperSpectrum(energy::EnergyScale, props::Dict{Symbol,Any}, arr::AxisArray; stagemap::Type{<:StageMapping}=DefaultStageMapping, livetime=get(props, :LiveTime, 1.0)  * ones(Float64, size(arr[2:end])...))
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        new{eltype(arr), ndims(arr) - 1, ndims(arr)}(
            arr,
            energy,
            livetime,
            props,
            stagemap,
            hash(arr, hash(energy, hash(livetime, hash(props, hash(stagemap)))))
        )
    end
    function HyperSpectrum(
        energy::EnergyScale,
        props::Dict{Symbol,Any},
        dims::Tuple,
        depth::Int,
        type::Type{<:Real};
        axisnames = ( "Y", "X", "Z", "A", "B", "C" ), 
        fov = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ),
        offset = ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ), #
        stagemap::Type{<:StageMapping}=DefaultStageMapping
    )
        axes = [ Axis{:Channel}(1:depth), #
            ( Axis{Symbol(axisnames[i])}(offset[i]:fov[i]/(dims[i]-1):offset[i]+fov[i]) for i in 1:length(dims) )... ]
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        counts = AxisArray(zeros(type, depth, dims...), axes...)
        livetime = fill(Float64(get(props,:LiveTime,1.0)), dims...)
        new{type, length(dims), length(dims)+1}(
            counts,
            energy,
            livetime,
            props,
            stagemap,
            hash(counts, hash(energy, hash(livetime, hash(props, hash(stagemap)))))
        )
    end
end

Base.hash(hs::HyperSpectrum, h::UInt) = hash(hs.hash, h)

function Base.show(io::IO, ::MIME"text/plain", hss::HyperSpectrum{T,N,NP}) where {T<:Real, N, NP}
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$T,$(axisnames(hss))}[$(name(hss)), $(hss.energy), $(depth(hss)) ch]",
    )
end

function Base.show(io::IO, ::MIME"text/html", hss::HyperSpectrum{T,N,NP}) where {T<:Real, N, NP}
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "<p>$sz HyperSpectrum{$T,$(axisnames(hss))}[$(name(hss))), $(hss.energy), $(depth(hss)) ch]</p>",
    )
end

function Base.show(io::IO, hss::HyperSpectrum{T,N,NP}) where {T<:Real, N, NP}
    sz = join(repr.(size(hss)), " × ")
    print(
        io,
        "$sz HyperSpectrum{$T,$(axisnames(hss))}[$(name(hss)), $(hss.energy), $(depth(hss)) ch]",
    )
end

Base.eltype(::HyperSpectrum{T,N,NP}) where {T<:Real,N,NP} = Spectrum{T}
Base.length(hss::HyperSpectrum) = prod(size(hss))
Base.ndims(::HyperSpectrum{T,N,NP}) where {T<:Real,N,NP} = N
Base.size(hss::HyperSpectrum) = size(hss.counts)[2:end]
Base.size(hss::HyperSpectrum, n) = size(counts(hss), n+1)
Base.axes(hss::HyperSpectrum) = Base.axes(counts(hss))[2:end]
Base.axes(hss::HyperSpectrum, n) = Base.axes(counts(hss), n+1)
Base.eachindex(hss::HyperSpectrum) = Base.OneTo(length(hss))
Base.stride(hss::HyperSpectrum, k::Int) = stride(counts(hss), k + 1)
Base.strides(hss::HyperSpectrum) = strides(counts(hss))[2:end]
depth(hss::HyperSpectrum) = size(counts(hss), 1)
NeXLCore.name(hss::HyperSpectrum) = get(hss.properties, :Name, "HyperSpectrum")
AxisArrays.axisnames(hss::HyperSpectrum) = AxisArrays.axisnames(hss.counts)[2:end]
AxisArrays.axisvalues(hss::HyperSpectrum) = AxisArrays.axisvalues(hss.counts)[2:end]
AxisArrays.axes(hss::HyperSpectrum) = AxisArrays.axes(hss.counts)[2:end]
AxisArrays.axes(hss::HyperSpectrum, dim::Int) = AxisArrays.axes(hss.counts,dim+1)
axisname(hss::HyperSpectrum, ax::Int) = AxisArrays.axisnames(hss.counts)[ax+1]
axisvalue(hss::HyperSpectrum, ax::Int, j::Int) = AxisArrays.axisvalues(hss.counts)[ax+1][j]
axisrange(hss::HyperSpectrum, ax::Int) = AxisArrays.axisvalues(hss.counts)[ax+1]
Base.CartesianIndices(hss::HyperSpectrum{T, N, NP}) where {T<:Real, N, NP} = CartesianIndices(size(hss))
Base.haskey(hss::HyperSpectrum, sym::Symbol) = haskey(hss.properties, sym)

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
    hss::HyperSpectrum{T,N,NP},
    ci::CartesianIndex
)::Spectrum{T} where { T<: Real, N, NP} 
    props = copy(hss.properties)
    props[:Cartesian] = ci
    props[:Name] = "$(name(hss))[$(join(["$i" for i in ci.I],","))]"
    props[:StageCoordinate]=coordinate(hss, ci)
    props[:LiveTime] = hss.livetime[ci]
    return Spectrum(hss.energy, counts(hss)[:, ci], props)
end
function Base.getindex(
    hss::HyperSpectrum{T, N, NP},
    idx::Vararg{Int,N},
)::Spectrum{T} where {T<:Real, N, NP}
    getindex(hss, CartesianIndex(idx...))
end
function Base.getindex(hss::HyperSpectrum{T, N, NP}, idx::Int)::Spectrum{T} where {T<:Real, N, NP}
    getindex(hss, CartesianIndices(hss)[idx])
end
function Base.getindex(
    hss::HyperSpectrum{T, N, NP},
    vidx::Vector
)::Vector{Spectrum{T}} where { T<: Real, N, NP} 
    return [ hss[i] for i in vidx ]
end
function Base.getindex(
    hss::HyperSpectrum{T, N, NP},
    rng::AbstractRange{<:Integer},
)::Vector{Spectrum{T}} where { T<: Real, N, NP} 
    return [ hss[i] for i in rng ]
end
function Base.getindex(
    hss::HyperSpectrum{T, N, NP},
    idx...
)::HyperSpectrum{T, N, NP} where {T<:Real, N, NP}
    props = copy(hss.properties)
    props[:Name] = "$(name(hss))[$idx]"
    return HyperSpectrum(hss.energy, props, counts(hss)[:, idx...], livetime=hss.livetime[idx...])
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
function Base.getindex(
    hss::HyperSpectrum,
    cxr1::CharXRay,
    cxr2::CharXRay
) 
    colorize(hss, [cxr1, cxr2], :All)
end
function Base.getindex(
    hss::HyperSpectrum,
    cxr1::CharXRay,
    cxr2::CharXRay,
    cxr3::CharXRay
) 
    colorize(hss, [cxr1, cxr2, cxr3], :Each)
end

function Base.get(hss::HyperSpectrum, sym::Symbol, def) 
    get(hss.properties, sym, def)
end


"""
    region(hss::HyperSpectrum{T, N}, ranges::AbstractRange...)::HyperSpectrum where {T<:Real,N}

Creates a view of a HyperSpectrum to represent the range of pixels within the `hss` HyperSpectrum.  Does
not copy the data or properties so any modifications to the region are also made to `hss`.
"""
function region(
    hss::HyperSpectrum{T, N, NP},
    ranges::AbstractRange...,
)::HyperSpectrum{T, N, NP} where {T<:Real, N, NP}
    HyperSpectrum(hss.energy, hss.properties, counts(hss)[:, ranges...], livetime=hss.livetime[ranges...])
end

function Base.getindex(hss::HyperSpectrum, sy::Symbol)
    @assert sy!=:LiveTime "Use the livetime(...) function to get the :LiveTime property for hyperspectra."
    getindex(hss.properties, sy)
end
function Base.setindex!(hss::HyperSpectrum, val::Any, sy::Symbol)
    @assert sy!=:LiveTime "Use the livetime!(...) function to set the :LiveTime property for hyperspectra."
    isequal(sy, :Detector) && @assert isequal(hss.energy, val.scale) "The energy calibration in $val must match the HyperSpectrum's: $(hss.energy)."
    setindex!(hss.properties, val, sy)
end

NeXLCore.energy(ch::Integer, hss::HyperSpectrum) = energy(ch, hss.energy)
channel(energy::Real, hss::HyperSpectrum) = channel(energy, hss.energy)
rangeofenergies(ch::Integer, hss::HyperSpectrum) =
    (energy(ch, hss.energy), energy(ch + 1, hss.energy))
properties(hss::HyperSpectrum)::Dict{Symbol,Any} = hss.properties

"""
	matching(spec::HyperSpectrum, resMnKa::Float64, lld::Int=1)::BasicEDS

Build an EDSDetector to match the channel count and energy scale in this spectrum.
"""
matching(spec::HyperSpectrum, resMnKa::Float64, lld::Int = 1)::BasicEDS =
    BasicEDS(depth(spec), spec.energy, MnKaResolution(resMnKa), lld)

matching(
    spec::HyperSpectrum,
    resMnKa::Float64,
    lld::Int,
    minByFam::Dict{Shell,Element},
)::BasicEDS = BasicEDS(depth(spec), spec.energy, MnKaResolution(resMnKa), lld, minByFam)


"""
    counts(hss::HyperSpectrum{T, N, NP})::Array{T,NP}

Creates type-friendly view of the counts data array.  Use of this function helps to avoid performance
penalties associated with boxing/unboxing the counts data.

    counts(hss::HyperSpectrum{T, N, NP}, ci::CartesianIndex)::Vector{T}

Access the counts data associated with the pixel `ci`.

    counts(hss::HyperSpectrum{T, N, NP}, ci::CartesianIndex, ch::Int)::T

Access the counts data at the pixel represented by `ci` and the channel represented by `ch`.
"""
counts(hss::HyperSpectrum{T, N, NP}) where {T<:Real, N, NP} = parent(hss.counts)
counts(hss::HyperSpectrum{T, N, NP}, ci::CartesianIndex, ch::Int) where {T<:Real, N, NP} = counts(hss)[ch, ci]
counts(hss::HyperSpectrum{T, N, NP}, ci::CartesianIndex) where {T<:Real, N, NP} = view(counts(hss)[:, ci])

"""
    compress(hss::HyperSpectrum)

Returns a HyperSpectrum with smaller or equal storage space to `hss` without losing or truncating any counts 
(note: AbstractFloat compresses to Float32 with loss of precision).  Can change the storage type and/or 
reduce the depth of hss.
"""
function compress(hss::HyperSpectrum{T, N, NP}) where {T<:Real, N, NP}
    data = counts(hss)
    (minval, maxval) = extrema(data)
    maxi = mapreduce(max, CartesianIndices(hss)) do ci
        something(findlast(c -> c != 0.0, data[:, ci]), 0)
    end
    data = maxi < depth(hss) ? data[1:maxi, axes(data)[2:end]...] : data
    if T isa Int16 || T isa Int32 || T isa Int64
        for newtype in filter(ty -> sizeof(T) > sizeof(ty), (Int8, Int16, Int32))
            if minval >= typemin(newtype) && maxval <= typemax(newtype)
                return HyperSpectrum(hss.energy, copy(hss.properties), newtype.(data), livetime=hss.livetime)
            end
        end
    elseif T isa UInt16 || T isa UInt32 || T isa UInt64
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

Sums a contiguous range of data channels into an Array. The dimension of the result is
one less than the dimension of the HyperSpectrum and is stored as a Float64 to ensure that not information is lost.
"""
function plane(
    hss::HyperSpectrum{T, N, NP},
    chs::AbstractUnitRange{<:Integer},
    normalize = false,
) where {T<:Real, N, NP}
    data = counts(hss)
    res = map(CartesianIndices(hss)) do ci
        sum(view(data, chs, ci), init=zero(Float64))
    end
    return normalize ? res ./= maximum(res) : res
end
function plane(
    hss::HyperSpectrum{T, N, NP},
    cxr::CharXRay,
    normalize = false,
) where {T<:Real, N, NP}
    plane(hss, fwhmroi(hss, cxr), normalize)
end


"""
   plane(hss::HyperSpectrum, ch::Int, normalize=false)

Extracts a single plane from a HyperSpectrum. The dimension of the result is
one less than the dimension of the HyperSpectrum.
"""
function plane(hss::HyperSpectrum, ch::Int, normalize = false)
    res = counts(hss)[ch, CartesianIndices(hss)]
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
    hss::HyperSpectrum{T, N, NP},
    cis::CartesianIndices,
    filt::Function = ci -> true,
)::Spectrum{T} where {T<:Real, N, NP}
    data, res = counts(hss), zeros(T, depth(hss))
    for ci in filter(ci -> filt(ci), cis), ch in eachindex(res)
        res[ch] = max(res[ch], data[ch, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    props[:Name] = "MaxPixel[$(props[:Name])]"
    return Spectrum(hss.energy, res, props)
end
function maxpixel(
    hss::HyperSpectrum{T, N, NP},
)::Spectrum{T} where {T<:Real, N, NP}
    data, res  = counts(hss), fill(zero(T), depth(hss))
    for ci in CartesianIndices(hss), ch in eachindex(res)
        res[ch] = max(res[ch], data[ch, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    props[:Name] = "MaxPixel[$(props[:Name])]"
    return Spectrum(hss.energy, res, props)
end

function minpixel(
    hss::HyperSpectrum{T, N, NP},
)::Spectrum{T} where {T<:Real, N, NP}
    data, res  = counts(hss), fill(typemax(T), depth(hss))
    for ci in CartesianIndices(hss), ch in eachindex(res)
        res[ch] = min(res[ch], data[ch, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    props[:Name] = "MaxPixel[$(props[:Name])]"
    return Spectrum(hss.energy, res, props)
end
minpixel(hss::HyperSpectrum, mask::BitArray) =
    minpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function minpixel(
    hss::HyperSpectrum{T, N, NP},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real, N, NP}
    data, res  = counts(hss), fill(typemax(T), depth(hss))
    for ci in filter(ci -> filt(ci), cis), ch in eachindex(res)
        res[ch] = min(res[ch], data[ch, ci])
    end
    props = copy(hss.properties)
    props[:LiveTime] = mean(hss.livetime)
    props[:Name] = "MaxPixel[$(props[:Name])]"
    return Spectrum(hss.energy, res, props)
end
avgpixel(hss::HyperSpectrum, mask::BitArray) =
    avgpixel(hss, CartesianIndices(hss), ci -> mask[ci])
function avgpixel(
    hss::HyperSpectrum{T, N, NP},
    cis::CartesianIndices,
    filt::Function = ci -> true,
) where {T<:Real, N, NP}
    data = counts(hss)
    res = zeros(T isa Int ? Int64 : Float64, depth(hss))
    for ci in filter(ci -> filt(ci), cis)
        res .+= @view data[:, ci]
    end
    res /= prod(size(hss))
    lt = mean(hss.livetime)
    props = copy(hss.properties)
    props[:Name] = something(name, "Average[$(props[:Name])]")
    return Spectrum(hss.energy, res, props)
end
function avgpixel(
    hss::HyperSpectrum{T, N, NP}
) where {T<:Real, N, NP}
    data = counts(hss)
    res = zeros(T isa Int ? Int64 : Float64, depth(hss))
    for ci in CartesianIndices(hss)
        res .+= @view data[:, ci]
    end
    res /= prod(size(hss))
    lt = mean(hss.livetime)
    props = copy(hss.properties)
    props[:Name] = something(name, "Average[$(props[:Name])]")
    return Spectrum(hss.energy, res, props)
end

"""
    sumcounts(hss::HyperSpectrum, cis::CartesianIndices = CartesianIndices(hss))

An array containing the number of counts at each pixel.
"""
function sumcounts(hss::HyperSpectrum, cis::CartesianIndices = CartesianIndices(hss))
    data = counts(hss)
    map(cis) do ci
        sum(view(data, :, ci), init=zero(Float64))
    end
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

function indexofmaxpixel(hss::HyperSpectrum{T, N, NP}, cis::CartesianIndices) where {T<:Real, N, NP}
    res, cix = zeros(T, depth(hss)), CartesianIndex[CartesianIndices(hss)[1] for _ = 1:depth(hss)]
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
    Base.sum(hss::HyperSpectrum{T, N}, func::Function) where {T<:Real, N}
    Base.sum(hss::HyperSpectrum{T, N}, cis::CartesianIndices) where {T<:Real, N}

Compute a sum spectrum for all or a subset of the pixels in `hss`.

Where func(hss::HyperSpectrum, ci::CartesianIndex)::Bool
"""
function Base.sum(
    hss::HyperSpectrum{T, N, NP},
    mask::Union{BitArray,Array{Bool}};
    name = nothing
) where {T<:Real, N, NP}
    @assert size(mask) == size(hss) "Mask size[$(size(mask))] ≠ Hyperspectrum size[$(size(hss))]"
    data = counts(hss)
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in CartesianIndices(hss)
        if mask[ci]
            res .+= data[:, ci]
            lt += hss.livetime[ci]
        end
    end
    props = copy(hss.properties)
    props[:Name] = something(name, "MaskedSum[$(props[:Name])]")
    return Spectrum(hss.energy, res, props)
end
function Base.sum(
    hss::HyperSpectrum{T, N, NP},
    filt::Function;
    name = nothing
) where {T<:Real, N, NP}
    data = counts(hss)
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in CartesianIndices(hss)
        if filt(hss,ci)
            res .+= data[:, ci]
            lt += hss.livetime[ci]
        end
    end
    props = copy(hss.properties)
    props[:Name] = something(name, "FilteredSum[$(props[:Name])]")
    return Spectrum(hss.energy, res, props)
end
function Base.sum(hss::HyperSpectrum{T, N, NP}; name = nothing) where {T<:Real, N, NP}
    data = counts(hss)
    res = zeros(T isa Int ? Int64 : Float64, depth(hss))
    for ci in CartesianIndices(hss)
        res .+= @view data[:, ci]
    end
    lt = sum(hss.livetime)
    props = copy(hss.properties)
    props[:Name] = something(name, "Sum[$(props[:Name])]")
    return Spectrum(hss.energy, res, props)
end
function Base.sum(hss::HyperSpectrum{T, N, NP}, cis::CartesianIndices; name = nothing) where {T<:Real, N, NP}
    data = counts(hss)
    res, lt = zeros(T isa Int ? Int64 : Float64, depth(hss)), 0.0
    for ci in cis
        res .+= data[:, ci]
        lt += hss.livetime[ci]
    end
    props = copy(hss.properties)
    props[:Name] = something(name, "IndexedSum[$(props[:Name])]")
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
    Gray.(N0f8.(plane(hss, chs, true)))

function fwhmroi(hss::HyperSpectrum, cxr::CharXRay)
    ecxr = energy(cxr)
    r = haskey(hss.properties, :Detector) ? resolution(ecxr, hss[:Detector]) : resolution(ecxr, MnKaResolution(130.0))
    return channel(ecxr-0.5*r, hss.energy):channel(ecxr+0.5*r, hss.energy)
end

"""
    roiimage(hss::HyperSpectrum, cxr::CharXRay)

Create a count map for the specified characteristic X-ray.  By default, integrates
for one FWHM at `cxr`.  If hss[:Detector] is an `EDSDetector`, the FWHM is taken from it.
Otherwise, a FWHM of 130 eV at Mn Kα is assumed.
"""
function roiimage(hss::HyperSpectrum, cxr::CharXRay)
    return roiimage(hss, fwhmroi(hss, cxr))
end

"""
    roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay})

Create an array of Gray images representing the intensity in each of the CharXRay lines
in `cxrs`.  They are normalized such the the most intense pixel in any of them defines white.
By default, integrates for one FWHM at `cxr`.  If hss[:Detector] is an `EDSDetector`, the FWHM 
is taken from it.   Otherwise, a FWHM of 130 eV at Mn Kα is assumed.
"""
function roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay})
    achs = map(cxr->fwhmroi(hss, cxr), cxrs)
    return roiimages(hss, achs)
end

"""
    colorize(krs::AbstractVector{KRatios}, red::Element, green::Element, blue::Element, normalize=:All[|:Each])

Create RGB colorized images from three `KRatios` or from three `Element`s.  The elements
are normalized relative to all `KRatios` in `krs`. The resulting images are scaled by the factor
`scale` to allow visualization of trace elements.
"""
function NeXLCore.colorize(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, normalize=:All)
    if normalize==:All
        imgs = roiimages(hss, cxrs)
        colorview(RGB, imgs[1], length(imgs) > 1 ? imgs[2] : zeroarray, length(imgs) > 2 ? imgs[3] : zeroarray )
    else
        # Normalize relative to max of each ROIs independently
        colorview(RGB, roiimage(hss,cxrs[1]), length(cxrs)>1 ? roiimage(hss,cxrs[2]) : zeroarray, length(cxrs) > 2 ? roiimage(hss, cxrs[3]) : zeroarray)
    end
end

"""
    block(hss::HyperSpectrum{T,N,NP}, steps::NTuple{N, Int})::HyperSpectrum{T,N,NP} where {T<:Real, N, NP}
    block(hss::HyperSpectrum{T,N,NP}, step::Int) where {T<:Real, N, NP}

Reduce the size of a HyperSpectrum by summing together blocks of adjacent pixels.  For example, `steps=(4,4)` would
sum together blocks of 16 spectra in `hss` to form a single pixel in the resulting `Hyperspectrum`.
"""
function block(hss::HyperSpectrum{T,N,NP}, step::Int) where {T<:Real, N, NP}
    block(hss, ntuple(_ -> step, N))
end
function block(hss::HyperSpectrum{T,N,NP}, steps::NTuple{N, Int})::HyperSpectrum where {T<:Real, N, NP}
    dims = size(hss) .÷ steps
    fov = [last(ax.val)-first(ax.val) for ax in AxisArrays.axes(hss.counts)[2:end]]
    offset = [first(ax.val) for ax in AxisArrays.axes(hss.counts)[2:end]]
    et = eltype(hss.counts) isa Integer ? widen(eltype(hss.counts)) : eltype(hss.counts)
    res = HyperSpectrum(hss.energy, copy(hss.properties), dims, depth(hss), et, 
        axisnames=axisnames(hss.counts)[2:end], fov=fov, offset=offset, stagemap=hss.stagemap)
    cx, rcx = counts(hss), counts(res) # Force the type to eliminate boxing
    ii = ntuple(_ -> 1, length(steps))
    tmp = zeros(et, depth(hss)) # accumulate it here
    for ci in CartesianIndices(dims) 
        tmp .= zero(et)
        lt = zero(eltype(res.livetime))
        for inner in CartesianIndices(steps)
            idx = (ci.I .- ii) .* steps .+ inner.I
            tmp .+= @view cx[:, idx...] # @view is important!!!
            lt += hss.livetime[idx...]
        end
        rcx[:, ci] .= tmp
        res.livetime[ci] = lt
    end
    res
end

"""
    linescan(hss::HyperSpectrum{T,2,3}, ci1::CartesianIndex{2}, ci2::CartesianIndex{2})

Extract pixels from `hss` along the line from `ci1` to `ci2` as a 1D HyperSpectrum.  Only works on 2D SpectrumImages.
"""
function linescan(hss::HyperSpectrum{T,2,3}, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}) where {T<:Real}
    fov = [ sqrt(sum(map(1:2) do ii
        a = AxisArrays.axes(hss, ii)
        (a[ci2.I[ii]]-a[ci1.I[ii]])^2
    end)) ]
    len = sum(abs.(ci2.I .- ci1.I)) + 1
    res = HyperSpectrum(hss.energy, copy(hss.properties), (len, ), depth(hss), T, 
        axisnames=[ "Linescan[$(ci1.I), $(ci2.I)]" ], fov = fov, offset = 0.0*fov)
    cxh, cxr, idx = counts(hss), counts(res), 0
    drawline(ci1, ci2, true) do ci
        idx+=1
        cxr[:, idx] .= @view cxh[:, ci...]
    end
    return res
end