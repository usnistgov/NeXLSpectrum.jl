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
        fov = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ], #
        offset = 0.0 * collect(fov), #
        stagemap::Type{<:StageMapping}=DefaultStageMapping, #
        livetime=fill(get(props, :LiveTime, 1.0), size(arr)[2:end]...)
    )
        axes = Axis[ Axis{:Channel}(1:size(arr,1)) ]
        for i in 1:ndims(arr)-1
            left = offset[i]-fov[i]/2.0
            push!(axes, Axis{Symbol(axisnames[i])}(left:fov[i]/(size(arr,i+1)-1):(left+fov[i]))) 
        end
        haskey(props, :Name) || (props[:Name] = "HS$(hyperspectrumCounter())")
        counts = AxisArray(arr, axes...)
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
        axes = Axis[ Axis{:Channel}(1:depth) ]
        for i in 1:length(dims)
            left = offset[i]-fov[i]/2.0
            push!(axes, Axis{Symbol(axisnames[i])}(left:fov[i]/(dims[i]-1):(left+fov[i]))) 
        end
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
Base.size(hss::HyperSpectrum) = size(hss.counts)[2:end]
Base.length(hss::HyperSpectrum) = prod(size(hss))
Base.ndims(::HyperSpectrum{T,N,NP}) where {T<:Real,N,NP} = N
Base.axes(hss::HyperSpectrum) = Base.axes(counts(hss))[2:end]
Base.firstindex(hss::HyperSpectrum) = 1
Base.lastindex(hss::HyperSpectrum) = length(hss)
Base.eachindex(hss::HyperSpectrum) = Base.OneTo(length(hss))
Base.strides(hss::HyperSpectrum) = strides(counts(hss))[2:end]

# Implement iterators over Spectrum{T}
Base.IteratorSize(::Type{<:HyperSpectrum{T, N, NP}}) where {T <: Real, N, NP} =
    Base.HasShape{NP}()
Base.IndexStyle(::Type{<:HyperSpectrum{T, N, NP}}) where {T <: Real, N, NP} = IndexCartesian()
Base.IteratorEltype(::Type{<:HyperSpectrum{T, N, NP}}) where {T <: Real, N, NP} = Base.HasEltype()
Base.iterate(hss::HyperSpectrum{T, N, NP}, state=1) where {T <: Real, N, NP} =
    state > length(hss) ? nothing : ( getindex(hss, state), state + 1 )
    
NeXLCore.name(hss::HyperSpectrum) = get(hss.properties, :Name, "HyperSpectrum")
AxisArrays.axisnames(hss::HyperSpectrum) = AxisArrays.axisnames(hss.counts)[2:end]
AxisArrays.axisvalues(hss::HyperSpectrum) = AxisArrays.axisvalues(hss.counts)[2:end]
AxisArrays.axes(hss::HyperSpectrum) = AxisArrays.axes(hss.counts)[2:end]
axisname(hss::HyperSpectrum, ax::Int) = AxisArrays.axisnames(hss.counts)[ax+1]
axisvalue(hss::HyperSpectrum, ax::Int, j::Int) = AxisArrays.axisvalues(hss.counts)[ax+1][j]
axisrange(hss::HyperSpectrum, ax::Int) = AxisArrays.axisvalues(hss.counts)[ax+1]
Base.CartesianIndices(hss::HyperSpectrum{T, N, NP}) where {T<:Real, N, NP} = CartesianIndices(size(hss))
Base.haskey(hss::HyperSpectrum, sym::Symbol) = haskey(hss.properties, sym)

depth(hss::HyperSpectrum) = size(counts(hss), 1)

function matches(spec1::HyperSpectrum, spec2::HyperSpectrum, tol::Float64 = 1.0)::Bool
    return abs(energy(1, spec1) - energy(1, spec2)) < tol * channelwidth(1, spec1) &&
           abs(energy(length(spec1), spec1) - energy(length(spec2), spec2)) <
           tol * channelwidth(length(spec1), spec1)
end

"""
    coordinate(hss::HyperSpectrum, idx::Tuple{<:Int})

Computes the stage coordinate centering the pixel specified by `idx` using the `StageMapping` `hss.stagemap`.
"""
function coordinate(hss::HyperSpectrum, idx::Tuple)
    av, an = axisvalues(hss.counts), axisnames(hss.counts)
    img_coord = Dict(  an[i+1] => av[i+1][ii] for (i,ii) in enumerate(idx) )
    stage_coord = get(hss.properties, :StagePosition, Dict(sy=>zero(typeof(av[2][1])) for sy in (:X, :Y) ))
    th = get(hss.properties, :ImageRotation, 0.0)
    return image2stage(hss.stagemap, stage_coord, img_coord, th)
end
coordinate(hss::HyperSpectrum, ci::CartesianIndex) = #
    coordinate(hss, ci.I)

"""
    dose(hss::HyperSpectrum) # Average dose per pixel
    dose(hss::HyperSpectrum, idx...) # Dose for the `idx` pixel

Returns the product of the probe current and the live-time on a per pixel basis.
"""
dose(hss::HyperSpectrum) = hss[:ProbeCurrent] * mean(hss.livetime)
dose(hss::HyperSpectrum, idx...)  = hss[:ProbeCurrent] * hss.livetime[idx...]

"""
    livetime!(hss::HyperSpectrum, lt::AbstractFloat, idx...)
    livetime!(hss::HyperSpectrum{T,N}, lt::AbstractFloat) # All pixels to lt

Set the livetime on a per pixel basis.
"""
function livetime!(hss::HyperSpectrum, lt::AbstractFloat, idx...) 
    hss.livetime[idx...] = lt
end
function livetime!(hss::HyperSpectrum, lt::AbstractFloat)  
    hss.livetime .= lt
end

livetime(hss::HyperSpectrum, idx...) = hss.livetime[idx...]

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
    fov = map(1:length(idx)) do i 
        ax = AxisArrays.axes(hss,i)
        ax[last(idx[i])]-ax[first(idx[i])]
    end
    offset = map(1:length(idx)) do i
        ax = AxisArrays.axes(hss,i)
        0.5*(ax[last(idx[i])]+ax[first(idx[i])])
    end 
    return HyperSpectrum(hss.energy, props, counts(hss)[:, idx...], livetime=hss.livetime[idx...], fov = fov, offset = offset)
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
channelwidth(ch::Int, hss::HyperSpectrum) = energy(ch + 1, hss) - energy(ch, hss)
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
reduce the depth of hss.  If there is nothing that can be done to reduce the size of `hss` then `compress(hss)===hss`.

Example:

     hs = compress(hs) # Replace `hs` with a version that uses less memory (if possible).
"""
function compress(hss::HyperSpectrum{T, N, NP}) where {T<:Integer, N, NP}
    (minval, maxval) = extrema(counts(hss))
    tyopts = minval < 0 ? (Int8, Int16, Int32, Int64) : (UInt8, UInt16, UInt32, UInt64)
    restype = tyopts[findfirst(tyopts) do ty
        minval >= typemin(ty) && maxval <= typemax(ty)
    end]
    @assert sizeof(restype) <= sizeof(T) "$restype"
    return compact(restype, hss)
end
function compress(hss::HyperSpectrum{T, N, NP}) where {T<:AbstractFloat, N, NP}
    (minval, maxval) = extrema(counts(hss))
    restype=T
    if sizeof(T) < sizeof(Float32) &&
        minval >= typemin(Float32) &&
        maxval <= typemax(Float32)
        restype=Float32
    end
    return compact(restype, hss)
end

"""
    compact(S::Type{<:Real}, hss::HyperSpectrum)

Converts the counts data type to S.  Also eliminates any channel planes on the end that
are entirely less than or equal to zero.
"""
function compact(S::Type{<:Real}, hss::HyperSpectrum)
    data = counts(hss)
    itr = Iterators.reverse(Base.axes(data,1))
    maxi = findfirst(itr) do ch # First non-negative channel
        !all((<=)(zero(eltype(data))), @view data[ch, Base.axes(data)[2:end]...])
    end
    @assert !isnothing(maxi) "There are no non-negative counts in this hyperspectrum."
    maxch = itr[maxi]
    if S != eltype(data) || maxch < size(data,1)
        sdata = maxch == size(data,1) ? S.(data) : S.(data[1:maxch, Base.axes(data)[2:end]...])
        return HyperSpectrum(hss.energy, copy(hss.properties), sdata, livetime=hss.livetime)
    else
        return hss
    end
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

Extracts a single channel plane from a HyperSpectrum. The dimension of the result is
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
(Accounts for differences in :LiveTime between pixels.)
"""
function roiimages(hss::HyperSpectrum, achs::AbstractVector{<:AbstractUnitRange{<:Integer}})
    data = counts(hss)
    ps = map(achs) do chs
        map(CartesianIndices(hss)) do ci
            sum(view(data, chs, ci), init=zero(Float64)) / hss.livetime[ci]
        end
    end
    maxval = maximum(map(p -> maximum(p), ps))
    return map(p -> Gray.(convert.(N0f8, p / maxval)), ps)
end

"""
    roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer})

Create a count map from the specified contiguous range of channels.
(Accounts for differences in :LiveTime between pixels.)
"""
function roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer})
    data = counts(hss)
    tmp = map(CartesianIndices(hss)) do ci
        sum(view(data, chs, ci), init=zero(Float64)) / hss.livetime[ci]
    end
    Gray.(N0f8.(tmp ./= maximum(tmp)))
end

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
roiimage(hss::HyperSpectrum, cxr::CharXRay) = roiimage(hss, fwhmroi(hss, cxr))

"""
    roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay})

Create an array of Gray images representing the intensity in each of the CharXRay lines
in `cxrs`.  They are normalized such the the most intense pixel in any of them defines white.
By default, integrates for one FWHM at `cxr`.  If hss[:Detector] is an `EDSDetector`, the FWHM 
is taken from it.   Otherwise, a FWHM of 130 eV at Mn Kα is assumed.
"""
roiimages(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}) = roiimages(hss, map(cxr->fwhmroi(hss, cxr), cxrs))

"""
    colorize(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, normalize=:All)

Create RGB colorized images from up to three `CharXRay` which define channels over which
the count signal is integrated. `normalize=:All` puts the intensities on a common scale
using the `roiimages(...)` method.  Otherwise each image is scaled independently based on
the brightest pixel.
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
    linescan(hss::HyperSpectrum{T,2,3}, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}, width::Int=1) 

Extract pixels from `hss` along the line from `ci1` to `ci2` as a 1D HyperSpectrum.  The `width` argument integrates
the linescan along a line perpendicular to the primary axis. Only works on 2D SpectrumImages and for odd values of `width`.  
The algorithm does double count the occasional pixel but the length of each perpendicular is maintained at `width`.
The `AxisArrays.axes(...)` scaling is maintained so lengths on the linescan can be compared to lengths on the original map.
"""
function linescan(hss::HyperSpectrum{T,2,3}, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}, width::Int=1) where {T<:Real}
    @assert width >= 1 "Only positive widths are supported in linescan."
    @assert ci1 in CartesianIndices(hss) "$ci1 is not inside $(CartesianIndices(hss))"
    @assert ci2 in CartesianIndices(hss) "$ci2 is not inside $(CartesianIndices(hss))"
    function perpray(c1, c2, w) # Computes a list of CI's starting at c1 towards c2 of length w skipping c1
        i, res = 0, CartesianIndex[]
        drawray(c1, c2) do ci
            if (i>0) && (i<=w)
                push!(res, CartesianIndex(ci))
            end
            (i+=1)<=w
        end
        return res
    end
    width % 2 == 0 && @info "Even widths are rounded down to the nearest odd width."
    fov = [ sqrt(sum(map(1:2) do ii
        a = AxisArrays.axes(hss, ii)
        (a[ci2.I[ii]]-a[ci1.I[ii]])^2
    end)) ]
    len = sum(abs.(ci2.I .- ci1.I)) + 1
    hs2 = HyperSpectrum(hss.energy, copy(hss.properties), (len, ), depth(hss), T, 
        axisnames=[ "Linescan[$(ci1.I), $(ci2.I), $width]" ], fov = fov, offset = 0.0*fov)
    cxh, cxr, idx = counts(hss), counts(hs2), 0
    drawline(ci1, ci2, true) do ci
        idx+=1
        cxr[:, idx] .= @view cxh[:, ci...]
        hs2.livetime[idx] = hss.livetime[ci...]
    end
    if (width-1) ÷ 2 > 0
        bounds = CartesianIndices(hss)
        for rot in ( [ 0 -1; 1 0 ], [ 0 1; -1 0 ])
            perp = rot * collect(ci2.I .- ci1.I)
            r1 = perpray(ci1, CartesianIndex((perp  .+ ci1.I)...), (width-1) ÷ 2)
            r2 = perpray(ci2, CartesianIndex((perp  .+ ci2.I)...), (width-1) ÷ 2)
            for (ci1r, ci2r) in zip(r1, r2)
                idx = 0
                drawline(ci1r, ci2r, true) do ci
                    idx+=1
                    if all(i->ci[i] in bounds.indices[i], 1:2)
                        cxr[:, idx] .+= @view cxh[:, ci...]
                        hs2.livetime[idx] += hss.livetime[ci...]
                    end
                end
            end
        end
    end
    return hs2
end

function gridize(img,step=16)
    for ci in CartesianIndices(img)
        if any(i->i%step==0, ci.I)
            img[ci] = 1-img[ci]
        end
    end
    return img 
end
