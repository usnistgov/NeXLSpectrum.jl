"""
Represents the processed spectral data necessary to efficiently filter-fit one or more unknown spectra.
A `FilterFitPacket` contains the data necessary to filter the unknown and to apply pre-filtered references.
If there are duplicate `FilteredReference` for an elemental ROI, the preference is for the first one.  This
allows you to fill in unavailable "FilteredReference" elemental ROIs with more general ones.
"""
struct FilterFitPacket{S<:Detector, T<:AbstractFloat}
    detector::S
    filter::TopHatFilter{T}
    references::Vector{FilteredReference{T}}

    function FilterFitPacket(det::S, filt::TopHatFilter{T}, refs::Vector{FilteredReference{T}})  where {S<:Detector, T<:AbstractFloat}
        # Only permit one FilteredReference for each type of label for a set of X-rays 
        filtrefs = FilteredReference{T}[]
        for ref in refs
            exists = false
            for fr in filtrefs
                # If a reference already exists of this type for these X-rays reassign it
                if (typeof(fr.label) == typeof(ref.label))  && (fr.label.xrays==ref.label.xrays)
                    # filtrefs[i] = ref
                    exists = true
                    break
                end
            end
            if !exists
                push!(filtrefs, ref)
            end
        end
        new{S,T}(det, filt, filtrefs)
    end
end # struct

function Base.show(io::IO, ffp::FilterFitPacket)
    xrays = join(map(r -> "\t$(r.label),", ffp.references), "\n")
    print(io, "References[\n\t$(ffp.detector), \n$xrays\n]")
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, ffp::FilterFitPacket)

Summarize the `FilteredReference` structs within a `FilterFitPacket` as a `DataFrame`.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, ffp::FilterFitPacket)
    function p2b(fref)
        croi = max(1,first(fref.roi)-first(fref.ffroi)):min(length(fref.data), last(fref.roi)-first(fref.ffroi))
        sum(fref.charonly) / (sum(fref.data[croi])-sum(fref.charonly))
    end 
    function s2n(fref)
        croi = max(1,first(fref.roi)-first(fref.ffroi)):min(length(fref.data), last(fref.roi)-first(fref.ffroi))
        sum(fref.charonly) / sqrt(sum(fref.data[croi])-sum(fref.charonly))
    end 
    DataFrame(
        :Lines => [ fr.label.xrays for fr in ffp.references],
        :Material => [ get(fr.label.spectrum, :Composition, nothing) for fr in ffp.references],
        :ROI => [ fr.roi for fr in ffp.references],
        # :FullROI => [ fr.ffroi for fr in ffp.references],
        Symbol("P-to-B") => p2b.(ffp.references), 
        Symbol("S-to-N") => s2n.(ffp.references) 
    )
end


"""
    missingReferences(ffp::FilterFitPacket, elms::Vector{Element}, e0::Float64, ampl=1.0e-5)

Returns a `Vector{Tuple{Vector{CharXRay}, UnitRange{Int64}}}` containing the ROIs for which a 
`FilteredReference` is missing from the `FilterFitPacket`.
"""
function missingReferences(ffp::FilterFitPacket, elms::Vector{Element}, e0::Float64, ampl=1.0e-5)
    # Find all collections of X-rays for which we'll need references
    allXrays = mapreduce(elm->NeXLSpectrum.labeledextents(elm, ffp.detector, ampl, e0), append!, elms)
    # Remove those for which there are references
    filter(allXrays) do xr
        isnothing(findfirst(ref->ref.label.xrays == xr[1], ffp.references))
    end
end

"""
Returns a tuple containing the unfiltered spectra associated with the references.
"""
spectra(ffp::FilterFitPacket) = (ref.label.spectrum for ref in ffp.references)

"""
A list of elements for which there are filtered references.
"""
NeXLCore.elms(ffp::FilterFitPacket) = union( [ element(ref.label) for ref in filter(r->r.label isa CharXRayLabel, ffp.references) ] )

struct ReferencePacket
    spectrum::Spectrum
    element::Element
    material::Material
end

Base.show(io::IO, rp::ReferencePacket) = print(io, "ReferencePacket[$(symbol(rp.element)), $(name(rp.material)), $(name(rp.spectrum))]")

"""
    reference( elm::Element, spec::Spectrum, mat::Material=spec[:Composition]; pc = nothing, lt = nothing, e0 = nothing, coating = nothing)::ReferencePacket
    reference(els::AbstractVector{Element}, spec::Spectrum, mat::Material = spec[:Composition]; pc = nothing, lt = nothing, e0 = nothing, coating = nothing)::Vector{ReferencePacket}
    
Construct a `ReferencePacket` from a `Spectrum` collected from the specified `Material` for the specified `Element`.
Often used with `references(...)` to build `FilterFitPacket`s.

Optional named arguments `pc`, `lt`, `e0`, `coating` allow you to specify the probe current, live time, beam energy and
sample coating.
"""
function reference(
    elm::Element,
    spec::Spectrum,
    mat::Material = spec[:Composition];
    pc = nothing,
    lt = nothing,
    e0 = nothing,
    coating = nothing,
)::ReferencePacket
    if !isnothing(lt)
        @assert lt > 0.0 "The live time must be larger than zero for $(spec[:Name])."
        spec[:LiveTime] = lt
    end
    if !isnothing(pc)
        @assert pc > 0.0 "The probe current must be larger than zero for $(spec[:Name])."
        spec[:ProbeCurrent] = pc
    end
    if !isnothing(e0)
        @assert e0 >= 1.0e3 "The beam energy ($(e0/1000.0)) keV is too low for $(spec[:Name])."
        spec[:BeamEnergy] = e0
    end
    if !isnothing(coating)
        @assert coating isa Film
        spec[:Coating] = coating
    end
    spec[:Composition] = mat
    @assert haskey(spec, :LiveTime) "The :LiveTime property must be defined for $(spec[:Name]).  (Use the `lt` keyword argument)"
    @assert haskey(spec, :ProbeCurrent) "The :ProbeCurrent property must be defined for $(spec[:Name]).  (Use the `pc` keyword argument)"
    @assert haskey(spec, :BeamEnergy) "The :BeamEnergy property must be defined for $(spec[:Name]).  (Use the `e0` keyword argument)"
    @assert haskey(mat, elm) "$(Symbol(elm)) is not present in $mat."
    return ReferencePacket(spec, elm, mat)
end

reference(elm::Element, filename::AbstractString, mat::Material; kwargs...) =
    reference(elm, loadspectrum(filename), mat; kwargs...)

reference(elm::Element, filename::AbstractString; kwargs...) =
    reference(elm, loadspectrum(filename); kwargs...)

function reference(
    els::AbstractVector{Element},
    spec::Spectrum,
    mat::Material = spec[:Composition];
    kwargs...)
    map(el->reference(el, spec, mat; kwargs...), els)
end
reference(elm::AbstractVector{Element}, filename::AbstractString, mat::Material; kwargs...) =
    reference(elm, loadspectrum(filename), mat; kwargs...)


"""
    references(refs::AbstractVector{ReferencePacket}, det::EDSDetector)::FilterFitPacket
    references(refs::AbstractVector{ReferencePacket}, fwhm::Float64)::FilterFitPacket

Constructs a FilterFitPacket from a vector of `ReferencePackets`.  Each `ReferencePacket` represents a 
single ROI for an element.  It is possible more than one `ReferencePacket` might be defined for an 
elemental ROI.  In this case, the `ReferencePacket` with the lower index will take preference over
later ones.  This allows you to fill in only the missing elemental ROIs using spectra collected from 
alternative materials.  For example, a spectrum from F₂Fe is suitable for the Fe K-lines but not the 
Fe L-lines. So we might specify F₂Fe first to specify the references for the Fe K-lines and then fill 
in the L-lines with a spectrum from pure Fe.
"""
function references(
    refs::AbstractVector{ReferencePacket},
    det::EDSDetector;
    ftype::Type{<:AbstractFloat} = Float64
)::FilterFitPacket
    chcount = det.channelcount
    @assert all(length(r.spectrum) == chcount for r in refs) "The number of spectrum channels must match the detector for all spectra."
    @assert length(refs) > 0 "Please provide at least one ReferencePacket in references(...)"
    # Build the top-hat filter for det
    ff = buildfilter(ftype, det)
    # Apply the top-hat filter to all refs. Trying to thread this fails. :-(
    frefs = mapreduce(append!, refs) do ref
        frefs = filterreference(ff, ref.spectrum, ref.element, ref.material)
        length(frefs)==0 && @warn "Unable to create any filtered ROI references for $(ref.element) from $(name(ref.material))."
        frefs
    end
    return FilterFitPacket(det, ff, frefs)
end
references(refs::AbstractVector{ReferencePacket}, fwhm::Float64; ftype::Type{<:AbstractFloat}=Float64) =
    references(refs, matching(first(refs).spectrum, fwhm), ftype=ftype)

fit_spectrum(spec::Spectrum, ffp::FilterFitPacket{S, T}) where { S<:Detector, T<: AbstractFloat } =
    fit_spectrum(FilteredUnknownW{T}, spec, ffp.filter, ffp.references)

fit_spectrum(specs::AbstractVector{<:Spectrum}, ffp::FilterFitPacket{S, T}) where { S<:Detector, T<: AbstractFloat } =
    ThreadsX.map(spec -> fit_spectrum(FilteredUnknownW{T}, spec, ffp.filter, ffp.references), specs)

"""
    fit_spectrum(spec::Spectrum, ffp::FilterFitPacket)::FilterFitResult
    fit_spectrum(specs::AbstractVector{<:Spectrum}, ffp::FilterFitPacket)::Vector{FilterFitResult}

Fit a `Spectrum` or a vector of `Spectrum` using the specified `FilterFitPacket`.  The result is a
`FilterFitResult` structure which contains k-ratios, residuals, etc. 


    fit_spectrum(hs::HyperSpectrum, ffp::FilterFitPacket; mode::Symbol=:Fast, zero = x -> max(0.0, x))::Array{KRatios}

  * `mode = :Fast` - Uses precomputed, filtered "vector" fit method.  No uncertainties are available.
  * `mode = :Intermediate` - Uses an optimized full fit without refits for negative k-ratios.
  * `mode = :Full` - Uses the full single spectrum fit algorithm including refitting when one or more k-ratio is less than zero.

Performs a filtered fit on a hyperspectrum returning an `Array{KRatios}`.

Selecting a mode:
  :Fast is good for generating k-ratio maps or exploratory analysis of a k-ratio map. :Full is best when a
  quantitative map of a high count hyperspectrum is desired.  Fit results for major elements are similar for
  all three but differ for minor and trace elements.  Particularly when a k-ratio is slightly negative. This
  negative k-ratio can effect the other k-ratios.  :Fast also works less well when many elements (>>10) (particularly
  interfering elements) are included in the fit. Unfortunately, :Intermediate and :Full slow down when many elements
  are fit - O(n(elements)²).

The following timing on a 512 x 512 x 2048 hyperspectrum fitting 21 elements with 33 distinct ROIs on a fast laptop
with 64 GiB memory give a relative feel for the speed of each algorithm.  Yes, :Fast is approximately 100x faster than
:Intermediate and used almost 80x less memory.  (64-bit references, 32-bit references use about 1/2 the memory and
take about 2/3 the time.)

| mode=         | Threads  |Run time (s)   | Allocations  | Memory (GiB) | GC time  |
|---------------|----------|---------------|--------------|--------------|----------|
| :Fast         | 4        | 13.8          | 2.11 M       | 4.72         | 4.2%     |
| :Intermediate | 4        | 1921.2        | 13.11 M      | 364.35       | 57.2%    |
| :Full         | 4        | 2641.6        | 6.16 G       | 860.6        | 54.7%    |
| :Fast         | 1        |           |  M       |          |      |
| :Intermediate | 1        |         |  M      |        |     |
| :Full         | 1        |         |  G       |         |     |

"""
function fit_spectrum(
    hs::HyperSpectrum,
    ffp::FilterFitPacket{S, T};
    mode::Symbol = :Fast,
    zero = x -> max(Base.zero(T), x),
    sigma = Base.zero(T)
)::Array{KRatios} where { S <: Detector, T <: AbstractFloat }
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector) "The detector for the hyper-spectrum must match the detector for the filtered references."
    @assert (mode==:Fast) || (mode==:Intermediate) || (mode==:Full) "The mode argument must be :Fast, :Intermediate or :Full."
    if mode == :Fast
        vq = VectorQuant(ffp.references, ffp.filter)
        return fit_spectrum(hs, vq, zero)
    elseif mode == :Intermediate
        return fit_spectrum_int(hs, ffp, zero)
    elseif mode == :Full
        return fit_spectrum_full(hs, ffp, sigma)
    end
    # Should never get here...
    return  Array{KRatio}[]
end

function fit_spectrum_int(
    hs::HyperSpectrum,
    ffp::FilterFitPacket{S, T},
    zero::Function
)::Array{KRatios} where { S <: Detector, T <: AbstractFloat }
    unklabel = UnknownLabel(hs)
    function _tophatfilterhs(unklabel, data, thf, scale) 
        @assert length(data) <= length(thf) "The reference spectra must have more channels than the hyperspectrum data."
        filtered = [filtereddatum(thf, data, i) for i in eachindex(data)]
        dp = [max(x, one(T)) for x in data] # To ensure covariance isn't zero or infinite precision
        covar = [filteredcovar(thf, dp, i, i) for i in eachindex(data)]
        FilteredUnknownW{T}(unklabel, scale, 1:length(data), data, filtered, covar)
    end
    fitcontiguousx(unk, ffs, chs) = #
        _buildscale(unk, ffs) * pinv(_buildmodel(ffs, chs), rtol = convert(T, 1.0e-6)) * extract(unk,chs)
    function _filterfitx(unk, ffs, fitrois) 
        Iterators.flatten(map(fitrois) do fr
            fitcontiguousx(
                unk,
                filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs),
                fr,
            )
        end)
    end
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector) "The detector for the hyper-spectrum must match the detector for the filtered references."
    krs = zeros(Float32, length(ffp.references), size(hs)...)
    len = 1:min(depth(hs),length(ffp.filter))
    data = zeros(T, length(ffp.filter))
    fitrois = ascontiguous(map(fd -> fd.ffroi, ffp.references))
    cdata = counts(hs)
    ThreadsX.foreach(CartesianIndices(hs)) do ci
        data[len] .= T.(cdata[len,ci])
        unk = _tophatfilterhs(unklabel, data, ffp.filter, convert(T,1/dose(hs,ci)))
        krs[:, ci] .= zero.(_filterfitx(unk, ffp.references, fitrois))
    end
    return ThreadsX.map(filter(ii -> ffp.references[ii].label isa CharXRayLabel, eachindex(ffp.references))) do i
        lbl = ffp.references[i].label
        rprops = properties(spectrum(lbl))
        KRatios(xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i, :, :])
    end
end

function fit_spectrum_full(
    hs::HyperSpectrum,
    ffp::FilterFitPacket{S, T},
    sigma::T
)::Array{KRatios} where { S <: Detector, T <: AbstractFloat }
    unklabel = UnknownLabel(hs)
    function _tophatfilterhs(unklabel, data, thf, scale) 
        @assert length(data) <= length(thf) "The reference spectra must have more channels than the hyperspectrum data."
        filtered = T[filtereddatum(thf, data, i) for i in eachindex(data)]
        dp = T[max(x, one(T)) for x in data] # To ensure covariance isn't zero or infinite precision
        covar = T[filteredcovar(thf, dp, i, i) for i in eachindex(data)]
        FilteredUnknownW{T}(unklabel, scale, 1:length(data), data, filtered, covar)
    end
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector) "The detector for the hyper-spectrum must match the detector for the filtered references."
    krs = zeros(Float32, length(ffp.references), size(hs)...)
    len = 1:depth(hs)
    data = zeros(T, length(ffp.filter))
    cdata = counts(hs)
    ThreadsX.foreach(CartesianIndices(hs)) do ci
        data[len] .= view(cdata, len, ci)
        unk = _tophatfilterhs(unklabel, data, ffp.filter, convert(T,1/dose(hs,ci)))
        uvs = _filterfit(unk, ffp.references, true)
        krs[:, ci] = map(ffp.references) do ref
            id = ref.label
            isnan(uvs, id) || (value(uvs, id) < sigma*σ(uvs, id)) ? Base.zero(eltype(krs)) : convert(eltype(krs), value(uvs, id))
        end
    end
    return ThreadsX.map(filter(ii -> ffp.references[ii].label isa CharXRayLabel, eachindex(ffp.references))) do i
        lbl = ffp.references[i].label
        rprops = properties(spectrum(lbl))
        KRatios( xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i, :, :])
    end
end