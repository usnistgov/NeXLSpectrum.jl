"""
Represents the processed spectral data necessary to efficiently filter-fit one or more unknown spectra.
A `FilterFitPacket` contains the data necessary to filter the unknown and to apply pre-filtered references.
If there are duplicate `FilteredReference` for an elemental ROI, the preference is for the first one.  This
allows you to fill in unavailable "FilteredReference" elemental ROIs with more general ones.
"""
struct FilterFitPacket{T<:Detector}
    detector::T
    filter::TopHatFilter
    references::Vector{FilteredReference}

    function FilterFitPacket(det::T, filt::TopHatFilter, refs::Vector{FilteredReference}) where {T <: Detector}
        # Only permit one FilteredReference for each type of label for a set of X-rays 
        filtrefs = FilteredReference[]
        for ref in refs
            exists = false
            for (i, fr) in enumerate(filtrefs)
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
        new{T}(det, filt, filtrefs)
    end
end # struct

function Base.show(io::IO, ffp::FilterFitPacket)
    lines = join(map(r -> "\t$(r.label),", ffp.references), "\n")
    print(io, "References[\n\t$(ffp.detector), \n$lines\n]")
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
NeXLCore.elms(ffp::FilterFitPacket) = union((elms(spec, true) for spec in spectra(ffp))...)

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
    det::EDSDetector,
)::FilterFitPacket
    chcount = det.channelcount
    @assert all(length(r.spectrum) == chcount for r in refs) "The number of spectrum channels must match the detector for all spectra."
    ff = buildfilter(det)
    
    frefs = if length(refs)>0
        mapreduce(append!, refs) do ref
            frefs = filterreference(ff, ref.spectrum, ref.element, ref.material)
            length(frefs)==0 && @warn "Unable to create any filtered ROI references for $(ref.element) from $(name(ref.material))."
            frefs
        end
    else
        FilteredReference[]
    end
    return FilterFitPacket(det, ff, frefs)
end
references(refs::AbstractVector{ReferencePacket}, fwhm::Float64)::FilterFitPacket =
    references(refs, matching(refs[1].spectrum, fwhm))

fit_spectrum(spec::Spectrum, ffp::FilterFitPacket) =
    fit_spectrum(FilteredUnknownW, spec, ffp.filter, ffp.references)

fit_spectrum(specs::AbstractVector{<:Spectrum}, ffp::FilterFitPacket) =
    map(spec -> fit_spectrum(FilteredUnknownW, spec, ffp.filter, ffp.references), specs)

"""
fit_spectrum(hs::HyperSpectrum, ffp::FilterFitPacket; mode::Symbol=:Fast, zero = x -> max(0.0, x))::Array{KRatios}

  * mode = :Fast - Uses precomputed, filtered "vector" fit method
  * mode = :Intermediate - Uses an optimized full fit without refits for k(s) < 0.0
  * mode = :Full - Uses the full single spectrum fit algorithm including refitting when one or more k < 0.0

Performs a filtered fit on a hyperspectrum returning an `Array{KRatios}`.

Selecting a mode:
  :Fast is good for generating k-ratio maps or exploratory analysis of a k-ratio map. :Full is best when a
  quantitative map of a high count hyperspectrum is desired.  Fit frefsults for major elements are similar for
  all three but differ for minor and trace elements.  Particularly when a k-ratio is slightly negative. This
  negative k-ratio can effect the other k-ratios.  :Fast also works less well when many elements (>>10) (particularly
  interfering elements) are included in the fit. Unfortunately, :Intermediate and :Full slow down when many elements
  are fit - O(n(elements)²).

The following timing on a 512 x 512 x 2048 hyperspectrum fitting 15 elements with 25 distinct ROIs on a fast laptop
with 64 GiB memory give a relative feel for the speed of each algorithm.  Yes, :Fast is approximately 20x faster than
:Intermediate and used almost 100x less memory.  (Single thread timings)

|---------------|----------------|--------------|--------------|----------|
| mode          | Run time (s)   | Allocations  | Memory (GiB) | GC time  |
|---------------|----------------|--------------|--------------|----------|
| :Fast         | 44.8           | 3.95 M       | 5.72         | 4.2%     |
| :Intermediate | 1305.1         | 24.83 M      | 523.7        | 2.7%     |
| :Full         | 2056           | 2.7 G        | 786.1        | 4.2%     |
|---------------|----------------|--------------|--------------|----------|
"""
function fit_spectrum(
    hs::HyperSpectrum,
    ffp::FilterFitPacket;
    mode::Symbol = :Fast,
    zero = x -> max(0.0, x),
    sigma = 0.0
)::Array{KRatios}
    _tophatfilterhs =
        (hs, data, thf, scale) -> begin
            @assert length(data) <= length(thf) "The reference spectra have fewer channels than the hyperspectrum data."
            filtered = Float64[filtereddatum(thf, data, i) for i in eachindex(data)]
            dp = Float64[max(x, 1.0) for x in data] # To ensure covariance isn't zero or infinite precision
            covar = Float64[filteredcovar(thf, dp, i, i) for i in eachindex(data)]
            return FilteredUnknownW(
                UnknownLabel(hs),
                scale,
                eachindex(data),
                data,
                filtered,
                covar,
            )
        end
    fitcontiguousx =
        (unk, ffs, chs) ->
            _buildscale(unk, ffs) *
            pinv(_buildmodel(ffs, chs), rtol = 1.0e-6) *
            extract(unk, chs)
    _filterfitx =
        (unk, ffs, fitrois) -> Iterators.flatten(
            map(
                fr -> fitcontiguousx(
                    unk,
                    filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs),
                    fr,
                ),
                fitrois,
            ),
        )
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector) "The detector for the hyper-spectrum must match the detector for the filtered references."
    if mode == :Fast
        vq = VectorQuant(ffp.references, ffp.filter)
        return fit_spectrum(vq, hs, zero)
    elseif mode == :Intermediate
        krs = zeros(Float32, length(ffp.references), size(hs)...)
        len = 1:depth(hs)
        data = zeros(Float64, length(ffp.filter))
        fitrois = ascontiguous(map(fd -> fd.ffroi, ffp.references))
        for ci in CartesianIndices(hs)
            data[len] = hs.counts[len, ci]
            unk = _tophatfilterhs(hs, data, ffp.filter, 1.0/dose(hs,ci))
            krs[:, ci] .= zero.(_filterfitx(unk, ffp.references, fitrois))
        end
        frefs = KRatios[]
        for i in filter(
            ii -> ffp.references[ii].label isa CharXRayLabel,
            eachindex(ffp.references),
        )
            k, lbl = krs[i], ffp.references[i].label
            rprops = properties(spectrum(lbl))
            push!(
                frefs,
                KRatios(
                    xrays(lbl),
                    properties(hs),
                    rprops,
                    rprops[:Composition],
                    krs[i, :, :],
                ),
            )
        end
        return frefs
    elseif mode == :Full
        krs = zeros(Float32, length(ffp.references), size(hs)...)
        len = 1:depth(hs)
        data = zeros(Float64, length(ffp.filter))
        for ci in CartesianIndices(hs)
            data[len] = hs.counts[len, ci]
            unk = _tophatfilterhs(hs, data, ffp.filter, 1.0/dose(hs,ci))
            uvs = _filterfit(unk, ffp.references, true)
            krs[:, ci] = map(ref.label for ref in ffp.references) do id
                if !(isnan(id, uvs) || (value(id, uvs) < sigma*σ(id, uvs)))
                    convert(Float32, value(id, uvs))
                else
                    0.0f0
                end
            end
        end
        frefs = KRatios[]
        for i in filter(
            ii -> ffp.references[ii].label isa CharXRayLabel,
            eachindex(ffp.references),
        )
            k, lbl = krs[i], ffp.references[i].label
            rprops = properties(spectrum(lbl))
            push!(
                frefs,
                KRatios(
                    xrays(lbl),
                    properties(hs),
                    rprops,
                    rprops[:Composition],
                    krs[i, :, :],
                ),
            )
        end
        return frefs
    else
        @assert false "The mode argument must be :Fast, :Intermediate or :Full."
    end
end


fit_spectrum(
    ty::Type{FilteredUnknownG},
    unk::Spectrum,
    ffp::FilterFitPacket,
    forcezeros::Bool = true,
) = fit_spectrum(ty, unk, ffp.filter, ffp.references, forcezeros)

fit_spectrum(
    ty::Type{FilteredUnknownG},
    unks::AbstractVector{<:Spectrum},
    ffp::FilterFitPacket,
    forcezeros::Bool = true,
) = fit_spectrum(ty, unks, ffp.filter, ffp.references, forcezeros)
