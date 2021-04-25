using Statistics
using DataAPI

abstract type FitResult end
"""
A measured set of correlated k-ratios.
"""
struct BasicFitResult <: FitResult
    label::Label
    kratios::UncertainValues
end

function findlabel(ffr::FitResult, cxr::CharXRay)
    lbls = labels(ffr.kratios)
    fa = findall(lbl -> (lbl isa CharXRayLabel) && (cxr in lbl.xrays), lbls)
    @assert length(fa) ≠ 0 "No label found for $cxr."
    @assert length(fa) == 1 "Multiple labels found for $cxr."
    return lbls[fa[1]]
end

Base.show(io::IO, ffr::FitResult) = print(io, "$(ffr.label)")
NeXLUncertainties.NeXLUncertainties.value(label::ReferenceLabel, ffr::FitResult) =
    NeXLUncertainties.value(label, ffr.kratios)
NeXLUncertainties.σ(label::ReferenceLabel, ffr::FitResult) = σ(label, ffr.kratios)
NeXLUncertainties.uncertainty(label::ReferenceLabel, ffr::FitResult) =
    uncertainty(label, ffr.kratios)
NeXLUncertainties.getindex(ffr::FitResult, label::ReferenceLabel) =
    getindex(ffr.kratios, label)
Base.keys(ffr::FitResult) = keys(ffr.kratios)
NeXLUncertainties.labels(ffr::FitResult) = labels(ffr.kratios)
function kratio(cxr::CharXRay, ffr::FitResult)
    lbl = findlabel(ffr, cxr)
    return uv(NeXLUncertainties.value(lbl, ffr), σ(lbl, ffr))
end
unknown(ffr::FitResult)::Label = ffr.label
Base.values(lbl::ReferenceLabel, ffrs::Vector{<:FitResult}) =
    [NeXLUncertainties.value(lbl, ffr.kratios, 0.0) for ffr in ffrs]
σs(lbl::ReferenceLabel, ffrs::Vector{<:FitResult}) =
    [σ(lbl, ffr.kratios, 0.0) for ffr in ffrs]

"""
    kratios(ffr::FitResult)::Vector{KRatio}

The k-ratios associated with each `CharXRayLabel` as a vector 'KRatio' objects.
"""
function kratios(ffr::FitResult)::Vector{KRatio}
    lbls = filter(labels(ffr.kratios)) do l
        (l isa CharXRayLabel) && haskey(properties(l)[:Composition], element(l))
    end
    return map(lbls) do lbl
        props = properties(lbl)
        KRatio(lbl.xrays, properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
    end
end

"""
    kratio(ffr::FitResult, cxr::CharXRay)::KRatio

Extract the k-ratio associated with the specified CharXRay (zero if not measured)
"""
function kratio(ffr::FitResult, cxr::CharXRay)
    lbls = labels(ffr.kratios)
    i = findfirst(lbl->(lbl isa CharXRayLabel) && (cxr in lbl.xrays), lbls)
    @assert !isnothing(i) "There is no k-ratio associated with $cxr."
    lbl, props = lbls[i], properties(lbls[i])
    KRatio(lbl.xrays, unknown(ffr).spectrum.properties, props, props[:Composition], ffr.kratios[lbl])
  end

"""
    extractStandards(ffr::FitResult, elm::Element, mat::Material)::Vector{KRatio}

Extract a `Vector{KRatio}` for `elm::Element` from a `ffr::FilterFitResult` measured from `mat::Material`. 
"""
function extractStandards(ffr::FitResult, elm::Element, mat::Material)::Vector{KRatio}
    @assert haskey(mat, elm) "The element $elm is not contained within $mat."
    lbls = filter(labels(ffr.kratios)) do l
        (l isa CharXRayLabel) && (elm == element(l))
    end
    krs = map(lbls) do lbl
        props = properties(lbl)
        kr = KRatio(lbl.xrays, properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
        kr.unkProps[:Composition] = mat
        kr
    end
    length(krs) == 0 && @info  "There are no k-ratios associated with the $elm in $ffr."
    return krs
end

"""
    extractStandard(ffr::FitResult, cxrs::AbstractVector{CharXRay}, mat::Material)::Vector{KRatio}

Extract a `KRatio` for the `CharXRay` from a `FilterFitResult` associated with
the `Material`.
"""
function extractStandards(ffr::FitResult, cxrs::AbstractVector{CharXRay}, mat::Material)::Vector{KRatio}
    lbls = filter(labels(ffr.kratios)) do lbl
        (lbl isa CharXRayLabel) && any(cxr->cxr in xrays(lbl), cxrs)
    end
    krs = map(lbls) do lbl
        props = properties(lbl)
        kr = KRatio(lbl.xrays, properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
        kr.unkProps[:Composition] = mat
        kr
    end
    length(krs) == 0 && @info "There are no k-ratios associated with any of $cxrs in $ffr."
    return krs
end

"""
    heterogeneity(lbl::ReferenceLabel, ffrs::Vector{FilterFitResult})

Computes the ratio of the standard deviation of the measured values over the mean calculated uncertainty
from the fit.  A value near 1 means the sample appears homogeneous and a value greater than 1 means the sample
appears heterogeneous.
"""
heterogeneity(lbl::ReferenceLabel, ffrs::Vector{<:FitResult}) =
    std(values(lbl, ffrs)) / mean(σs(lbl, ffrs))

"""
    NeXLUncertainties.asa(::Type{DataFrame}, ffrs::AbstractVector{<:FitResult}; charOnly = true, withUnc = false, pivot = false)

Return generic `FitResult` as a DataFrame.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    ffrs::AbstractVector{<:FitResult};
    charOnly = true,
    withUnc = false,
    pivot = false,
)::DataFrame
    lbls = sort(collect(Set(mapreduce(labels, append!, ffrs))))
    lbls = charOnly ? filter(lbl -> lbl isa CharXRayLabel, lbls) : lbls
    if pivot
        res = DataFrame(ROI = lbls)
        for ffr in ffrs
            vals = [NeXLUncertainties.value(lbl, ffr.kratios, missing) for lbl in lbls]
            insertcols!(res, ncol(res) + 1, Symbol(repr(ffr.label)) => vals)
            if withUnc
                unc = [σ(lbl, ffr.kratios, missing) for lbl in lbls]
                insertcols!(res, ncol(res) + 1, Symbol('Δ' * repr(ffr.label)) => unc)
            end
        end
    else
        rowLbls = [repr(ffr.label) for ffr in ffrs]
        res = DataFrame(Symbol("Spectra") => rowLbls)
        for lbl in lbls
            vals = [NeXLUncertainties.value(lbl, ffr.kratios, missing) for ffr in ffrs]
            insertcols!(res, ncol(res) + 1, Symbol(repr(lbl)) => vals)
            if withUnc
                unc = [σ(lbl, ffr.kratios, missing) for ffr in ffrs]
                insertcols!(res, ncol(res) + 1, Symbol('Δ' * repr(lbl)) => unc)
            end
        end
    end
    return res
end

function NeXLUncertainties.asa(::Type{DataFrame}, fr::FitResult; withUnc = false)
    withUnc ?
    DataFrame(
        ROI = labels(fr),
        k = map(l -> NeXLUncertainties.value(l, fr.kratios), labels(fr)),
        δk = map(l -> σ(l, fr.kratios), labels(fr)),
    ) : #
    DataFrame(
        ROI = labels(fr),
        k = map(l -> NeXLUncertainties.value(l, fr.kratios), labels(fr)),
    )
end

function DataAPI.describe(ffrs::Vector{<:FitResult})::DataFrame
    df = asa(DataFrame, ffrs)[:, 2:end]
    desc = describe(df, :mean, :std, :min, :q25, :median, :q75, :max)
    lbls = filter(
        lbl -> lbl isa CharXRayLabel,
        sort(collect(Set(mapreduce(labels, append!, ffrs)))),
    )
    insertcols!(desc, 4, :hetero => [heterogeneity(lbl, ffrs) for lbl in lbls])
    return desc
end

"""
    FilterFitResult

Represents the result of fitting either FilteredUnknownG or FilteredUnknownW to a FilteredUnknown.

Struct elements

    label::UnknownLabel  # Identifies the unknown
    kratios::UncertainValues # Labeled with ReferenceLabel objects
    roi::UnitRange{Int} # Range of channels fit
    raw::Vector{Float64} # Raw spectrum data
    residual::Vector{Float64} # Residual spectrum
"""
struct FilterFitResult <: FitResult
    label::UnknownLabel
    kratios::UncertainValues
    roi::UnitRange{Int}
    raw::Vector{Float64}
    residual::Vector{Float64}
    peakback::Dict{ReferenceLabel,NTuple{2,Float64}}
end

"""
    residual(ffr::FilterFitResult)::Spectrum

A Spectrum containing the histogram representing the unknown spectrum
minus the fitted characteristic peaks shapes times the best fit coefficient.
"""
function residual(ffr::FilterFitResult)::Spectrum
    props = copy(properties(ffr.label.spectrum))
    props[:Name] = "Residual[$(props[:Name])]"
    return Spectrum(ffr.label.spectrum.energy, ffr.residual, props)
end

spectrum(ffr::FilterFitResult)::Spectrum = ffr.label.spectrum

"""
    characteristiccounts(ffr::FiterFitResult, strip)

Number of spectrum counts that were accounted for by the fitted elements with the `strip` Element(s) removed.
"""
characteristiccounts(ffr::FilterFitResult, strip) =
    sum(element(ref) in strip ? 0.0 : v[1] - v[2] for (ref, v) in ffr.peakback) # sum(ffr.raw[ffr.roi]-ffr.residual[ffr.roi])

"""
    peaktobackground(ffr::FilterFitResult, backwidth::Float64=10.0)::Float64

The peak-to-background ratio as determined from the raw and residual spectra integrated over
the fit region-of-interest and scaled to `backwidth` eV of continuum (nominally 10 eV).
"""
function peaktobackground(
    ffr::FilterFitResult,
    klabel::ReferenceLabel,
    backwidth::Float64 = 10.0,
)::Float64
    unk = spectrum(unknown(ffr))
    peak, back = ffr.peakback[klabel]
    return (peak * (energy(klabel.roi.stop, unk) - energy(klabel.roi.start, unk))) /
           (backwidth * back)
end


"""
    NeXLUncertainties.asa(::Type{DataFrame}, ffr::FilterFitResult; charOnly::Bool=true, material=nothing)::DataFrame

Tabulate details about each region-of-interest in the 'FilterFitResult' in a 'DataFrame'.
  * If charOnly then only display characteristic X-ray data (not escapes etc.)
  * If `material` is a Material then the computed k-ratio will also be tabulated along with kmeas/kcalc.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    ffr::FilterFitResult;
    charOnly::Bool = true,
    material = nothing,
)::DataFrame
    lbls, klbl, std, kr, dkr, roi1, roi2, peak, back, ptob = UnknownLabel[],
    ReferenceLabel[],
    String[],
    Float64[],
    Float64[],
    Int[],
    Int[],
    Float64[],
    Float64[],
    Float64[]
    kcalc, ratio = Float64[], Float64[]
    for lbl in NeXLUncertainties.sortedlabels(ffr.kratios)
        if (!charOnly) || (lbl isa CharXRayLabel)
            stdspec, unkspec = properties(lbl), ffr.label.spectrum
            push!(lbls, ffr.label)
            push!(std, stdspec[:Name])
            push!(klbl, lbl)
            push!(roi1, lbl.roi.start)
            push!(roi2, lbl.roi.stop)
            push!(kr, NeXLUncertainties.value(lbl, ffr.kratios))
            push!(dkr, σ(lbl, ffr.kratios))
            pb = ffr.peakback[lbl]
            push!(peak, pb[1])
            push!(back, pb[2])
            push!(ptob, peaktobackground(ffr, lbl))
            if material isa Material
                kc = NaN64
                if lbl isa CharXRayLabel
                    try
                        zcu = zafcorrection(
                            XPP,
                            ReedFluorescence,
                            NullCoating,
                            material,
                            lbl.xrays,
                            unkspec[:BeamEnergy],
                        )
                        zcs = zafcorrection(
                            XPP,
                            ReedFluorescence,
                            NullCoating,
                            stdspec[:Composition],
                            lbl.xrays,
                            stdspec[:BeamEnergy],
                        )
                        kc = k(zcu, zcs, unkspec[:TakeOffAngle], stdspec[:TakeOffAngle])
                    catch
                        kc = NaN64
                    end
                end
                push!(kcalc, kc)
                push!(ratio, NeXLUncertainties.value(lbl, ffr.kratios) / kc)
            end
        end
    end
    return material isa Material ?
           DataFrame(
        Spectrum = lbls,
        Feature = klbl,
        Reference = std,
        Start = roi1,
        Stop = roi2,
        K = kr,
        dK = dkr,
        Peak = peak,
        Back = back,
        PtoB = ptob,
        Kcalc = kcalc,
        Ratio = ratio,
    ) :
           DataFrame(
        Spectrum = lbls,
        Feature = klbl,
        Reference = std,
        Start = roi1,
        Stop = roi2,
        K = kr,
        dK = dkr,
        Peak = peak,
        Back = back,
        PtoB = ptob,
    )
end

"""
    filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::AbstractVector{FilteredReference})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
function filteredresidual(
    fit::FilterFitResult,
    unk::FilteredUnknown,
    ffs::AbstractVector{FilteredReference},
)::Vector{Float64}
    scaled(ff) =
        (NeXLUncertainties.value(ff.label, fit) * (ff.scale / unk.scale)) *
        extract(ff, unk.ffroi)
    return unk.filtered - mapreduce(scaled, +, ffs)
end

