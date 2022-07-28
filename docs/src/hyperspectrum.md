# ![](NeXL_sm.png)Spectrum
## [Working with HyperSpectrum objects](@id hyperspectrum_methods)

`HyperSpectrum{T<:Real, N, NP} <: AbstractArray{Spectrum{T}, N}` represents an N-dimensional array of `Spectrum{T}`
which share common properties. N = 1 represents a line scan, N = 2 represents a spectrum image, N = 3 can represent
a slice-and-view type hyperspectrum cube or a time series of spectrum images.  Higher N are also possible.

To construct an empty `HyperSpectrum` use
```julia
hs = HyperSpectrum(
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
```

To construct a `HyperSpectrum` from an existing array of (N+1) dimensions.  The first dimension is the channel data like (ch, Y, X, ...)
```julia
hs = HyperSpectrum(
    energy::EnergyScale, 
    props::Dict{Symbol,Any}, 
    arr::Array{<:Real};
    axisnames = ( "Y", "X", "Z", "A", "B", "C" ), #
    fov = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ], #
    offset = zeros(length(fov)), #
    stagemap::Type{<:StageMapping}=DefaultStageMapping, #
    livetime=fill(get(props, :LiveTime, 1.0), size(arr)[2:end]...)
)
```

We will assume a 2D hyperspectrum ("spectrum image") for these examples.

The `Spectrum` representing individual pixels in a `HyperSpectrum` are accessed using
```julia
hs[10,20]  # access the 10th row, 20th column
```
To index energy planes within the `HyperSpectrum` you can use
```julia
plane(
    hss::HyperSpectrum{T, N, NP},
    chs::AbstractUnitRange{<:Integer},
    normalize = false,
)
```
which sums over `chs` and optionally normalizes to the brightest pixel.
Alternatively,
```julia
plane(
    hss::HyperSpectrum{T, N, NP},
    cxr::CharXRay,
    normalize = false,
)
```
will extract one FWHM centered on `cxr`.  Like `plane(hss,n"Fe K-L3")`.

To extract `HyperSpectrum` data items representing regions from the `HyperSpectrum`
```julia
hs[10:20,20:30] # extract a HyperSpectrum representing the 10th to 20th row and 20th to 30th columns.
hs[10:2:20,20:2:30] # extract a HyperSpectrum representing every other pixel in the 10th to 20th row and 20th to 30th columns.
```

To extract a linescan from within a 2D `HyperSpectrum`
```julia
linescan(hss::HyperSpectrum{T,2,3}, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}, width::Int=1)
```

Mostly, a `HyperSpectrum` acts like an `Array{Spectrum}` and you can view and process the 
individual spectra this way.  However, this is often inefficient as a new `Spectrum` datum
must be allocated each coordinate.  There are functions designed to operate more efficiently
on the channel data.

```julia
sum(hs) # Produce a sum spectrum
maxpixel(hs) # produce a Bright's max-pixel spectrum
```
To find out which pixel contains the max-pixel for each channel in the `HyperSpectrum`
```julia
indexofmaxpixel(hss::HyperSpectrum) # An array of CartesianIndices
indexofmaxpixel(hss::HyperSpectrum, ch::Int) # A CartesianIndices for channel `ch`
```

To visualize planes within the `HyperSpectrum` as images
```julia
hss[n"Fe K-L3"]  # A FWHM ROI around the Fe K-L3 transition
roiimage(hss::HyperSpectrum, chs::AbstractUnitRange{<:Integer})
```
For an RGB image of up to three lines
```julia
colorize(hss::HyperSpectrum, cxrs::AbstractVector{CharXRay}, normalize=:All)
colorize(hs, [ n"Fe K-L3", n"Si K-L3", n"Al K-L3" ])
```

Quantifying `HyperSpectra` is like quantifying spectra.   There is a specialization 
of `fit_spectrum(...)` and `quantify(...)` optimized for hyper-spectra.

```julia
fit_spectra(
    hs::HyperSpectrum,
    ffp::FilterFitPacket{S, T};
    mode::Symbol = :Fast[|:Intermediate|:Full]
    zero = x -> max(Base.zero(T), x),
    sigma = Base.zero(T)
)::Array{KRatios}
```
`:Fast` uses a highly optimized but less precise algorithm that works for up to 
approximately 15 to 20 elements.  `:Intermediate` and `:Full` perform the slower
weighted filtered least-square fit.  `:Intermediate` sets k-ratios less than zero
to zero but doesn't refit, while `:Full` removes k-ratios that are less than zero 
and refits.

To convert the `Array{KRatios}` from `fit_spectra` to `Array{Material}`
```julia
NeXLMatrixCorrection.quantify(
    measured::AbstractVector{<:KRatios};
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = NullFluorescence,
    cc::Type{<:CoatingCorrection} = NullCoating,
    kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
    coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
)
```