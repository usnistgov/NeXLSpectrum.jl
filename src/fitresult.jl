using Statistics
using DataAPI
using Procrastinate

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

Base.show(io::IO, ffr::FitResult) = print(io, "FitResult($(ffr.label))")
NeXLUncertainties.NeXLUncertainties.value(ffr::FitResult, label::ReferenceLabel) =
    NeXLUncertainties.value(ffr.kratios, label)
NeXLUncertainties.σ(ffr::FitResult, label::ReferenceLabel) = σ(ffr.kratios, label)
NeXLUncertainties.uncertainty(ffr::FitResult, label::ReferenceLabel) =
    uncertainty(ffr.kratios, label)
NeXLUncertainties.getindex(ffr::FitResult, label::ReferenceLabel) =
    getindex(ffr.kratios, label)
Base.keys(ffr::FitResult) = keys(ffr.kratios)
NeXLUncertainties.labels(ffr::FitResult) = labels(ffr.kratios)
function kratio(cxr::CharXRay, ffr::FitResult)
    lbl = findlabel(ffr, cxr)
    return uv(NeXLUncertainties.value(ffr, lbl), σ(ffr, lbl))
end
unknown(ffr::FitResult)::Label = ffr.label
Base.values(ffrs::Vector{<:FitResult}, lbl::ReferenceLabel) =
    [NeXLUncertainties.value(ffr.kratios, lbl, 0.0) for ffr in ffrs]
σs(ffrs::Vector{<:FitResult}, lbl::ReferenceLabel) =
    [σ(ffr.kratios, lbl, 0.0) for ffr in ffrs]

"""
    kratios(ffr::FitResult)::Vector{KRatio}

The k-ratios associated with each `CharXRayLabel` as a vector 'KRatio' objects.
"""
function kratios(ffr::FitResult)::Vector{KRatio}
    lbls = filter(labels(ffr.kratios)) do l
        (l isa CharXRayLabel) && haskey(NeXLCore.properties(l)[:Composition], element(l))
    end
    return map(lbls) do lbl
        props = NeXLCore.properties(lbl)
        KRatio(lbl.xrays, NeXLCore.properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
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
    lbl, props = lbls[i], NeXLCore.properties(lbls[i])
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
        props = NeXLCore.properties(lbl)
        kr = KRatio(lbl.xrays, NeXLCore.properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
        kr.unkProps[:Composition] = mat
        kr
    end
    length(krs) == 0 && @info  "There are no k-ratios associated with the $elm in $ffr."
    return krs
end

"""
    extractStandards(ffr::FitResult, cxrs::AbstractVector{CharXRay}, mat::Material)::Vector{KRatio}

Extract a `KRatio` for the `CharXRay` from a `FilterFitResult` associated with
the `Material`.
"""
function extractStandards(ffr::FitResult, cxrs::AbstractVector{CharXRay}, mat::Material)::Vector{KRatio}
    lbls = filter(labels(ffr.kratios)) do lbl
        (lbl isa CharXRayLabel) && any(cxr->cxr in xrays(lbl), cxrs)
    end
    krs = map(lbls) do lbl
        props = NeXLCore.properties(lbl)
        kr = KRatio(lbl.xrays, NeXLCore.properties(unknown(ffr)), props, props[:Composition], ffr.kratios[lbl])
        kr.unkProps[:Composition] = mat
        kr
    end
    length(krs) == 0 && @info "There are no k-ratios associated with any of $cxrs in $ffr."
    return krs
end

"""
    heterogeneity(ffrs::Vector{FilterFitResult}, lbl::ReferenceLabel)

Computes the ratio of the standard deviation of the measured values over the mean calculated uncertainty
from the fit.  A value near 1 means the sample appears homogeneous and a value greater than 1 means the sample
appears heterogeneous.
"""
heterogeneity(ffrs::Vector{<:FitResult}, lbl::ReferenceLabel) =
    std(values(ffrs, lbl)) / mean(σs(ffrs, lbl))


"""
    FilterFitResult

Represents the result of fitting a FilteredUnknownW to a FilteredUnknown.

Struct elements

    label::UnknownLabel  # Identifies the unknown
    kratios::UncertainValues # Labeled with ReferenceLabel objects
    roi::UnitRange{Int} # Range of channels fit
    raw::Vector{T} # Raw spectrum data
    residual::Deferred # Residual spectrum
    peakback::Deferred # {Dict{ReferenceLabel,NTuple{3,Float64}}} with peak counts, background counts and counts/(nAs)


Use `import DataFrames; DataFrame(ffr::FilterFitResult)` to summarize in tabular form.
"""
struct FilterFitResult{T <: AbstractFloat} <: FitResult
    label::UnknownLabel
    kratios::UncertainValues
    roi::UnitRange{Int}
    raw::Vector{T}
    residual::Deferred
    peakback::Deferred
end

Base.show(io::IO, ffr::FilterFitResult) = print(io, "FitResult($(ffr.label))")


"""
    residual(ffr::FilterFitResult)::Spectrum

A Spectrum containing the histogram representing the unknown spectrum
minus the fitted characteristic peaks shapes times the best fit coefficient.
"""
residual(ffr::FilterFitResult) = ffr.residual()

"""
    spectrum(ffr::FilterFitResult)::Spectrum

Returns the original unknown spectrum.
"""
spectrum(ffr::FilterFitResult) = ffr.label.spectrum

"""
    characteristiccounts(ffr::FiterFitResult, strip)

Number of spectrum counts that were accounted for by the fitted elements with the `strip` Element(s) removed.
"""
characteristiccounts(ffr::FilterFitResult, strip) =
    sum(element(ref) in strip ? 0.0 : v[1] - v[2] for (ref, v) in ffr.peakback())

"""
    peaktobackground(ffr::FilterFitResult, backwidth::Float64=10.0)::Float64

The peak-to-background ratio as determined from the raw and residual spectra integrated over
the fit region-of-interest and scaled to `backwidth` eV of continuum (nominally 10 eV).
"""
function peaktobackground(
    ffr::FilterFitResult{T},
    klabel::ReferenceLabel,
    backwidth::T = convert(T, 10.0),
) where { T <: AbstractFloat }
    unk, pb = spectrum(unknown(ffr)), ffr.peakback()
    peak, back, _ = pb[klabel]
    return (peak * (energy(klabel.roi.stop, unk) - energy(klabel.roi.start, unk))) /
           (backwidth * back)
end


"""
    filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::AbstractVector{FilteredReference})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
function filteredresidual(
    fit::FilterFitResult{T},
    unk::FilteredUnknown,
    ffs::AbstractVector{FilteredReference{T}},
) where { T <: AbstractFloat }
    return unk.filtered - mapreduce(+, ffs) do ff
        (value(fit, ff.label) * (ff.scale / unk.scale)) * extract(ff, unk.ffroi)
    end
end


"""
    NeXLMatrixCorrection.estimatecoating(fr::FitResult, substrate::Material, coating::Material, cxr::CharXRay, mc::Type{<:MatrixCorrection}=XPP)::Film

Estimate the mass-thickness of `coating` on bulk `substrate` using the X-ray `cxr` for which there is KRatio in `fr`.  The result is
assigned to the properties of `fr` so that it can be account for in the matrix correction process.

This method is intended for use on standards where the substrate composition is known a priori.
"""
function NeXLMatrixCorrection.estimatecoating(fr::FitResult, substrate::Material, coating::Material, cxr::CharXRay, mc::Type{<:MatrixCorrection}=XPP)::Film
    NeXLCore.properties(fr.label.spectrum)[:Coating] = estimatecoating(substrate, coating, kratio(fr, cxr), mc)
end