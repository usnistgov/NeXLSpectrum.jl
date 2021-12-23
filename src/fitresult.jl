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
    extractStandards(ffr::FitResult, cxrs::AbstractVector{CharXRay}, mat::Material)::Vector{KRatio}

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
    heterogeneity(ffrs::Vector{FilterFitResult}, lbl::ReferenceLabel)

Computes the ratio of the standard deviation of the measured values over the mean calculated uncertainty
from the fit.  A value near 1 means the sample appears homogeneous and a value greater than 1 means the sample
appears heterogeneous.
"""
heterogeneity(ffrs::Vector{<:FitResult}, lbl::ReferenceLabel) =
    std(values(ffrs, lbl)) / mean(σs(ffrs, lbl))

"""
    NeXLUncertainties.asa(::Type{DataFrame}, ffrs::AbstractVector{<:FitResult}; charOnly = true, withUnc = false, pivot = false)

Return generic `FitResult` as a DataFrame.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    ffrs::AbstractVector{<:FitResult};
    charOnly = true,
    withUnc = false,
    format = :normal # :pivot or :long
)::DataFrame
    lbls = sort(collect(Set(mapreduce(labels, append!, ffrs))))
    lbls = charOnly ? filter(lbl -> lbl isa CharXRayLabel, lbls) : lbls
    if format==:pivot
        res = DataFrame(ROI = lbls)
        for ffr in ffrs
            vals = [value(ffr.kratios, lbl, missing) for lbl in lbls]
            insertcols!(res, ncol(res) + 1, Symbol(repr(ffr.label)) => vals)
            if withUnc
                unc = [σ(ffr.kratios, lbl, missing) for lbl in lbls]
                insertcols!(res, ncol(res) + 1, Symbol('Δ' * repr(ffr.label)) => unc)
            end
        end
    elseif format==:long
        if withUnc
            res = DataFrame(Spectrum=String[], ROI = String[], k = Float64[], dk = Float64[]) 
            for ffr in ffrs, lbl in lbls
                push!(res, ( repr(ffr.label), repr(lbl), value(ffr.kratios, lbl, missing), σ(ffr.kratios, lbl, missing) ))
            end
        else
            res = DataFrame(Spectrum=String[], ROI = String[], k = Float64[]) 
            for ffr in ffrs, lbl in lbls
                push!(res, ( repr(ffr.label), repr(lbl), value(ffr.kratios, lbl, missing)))
            end
        end
    else
        rowLbls = [repr(ffr.label) for ffr in ffrs]
        res = DataFrame(Symbol("Spectra") => rowLbls)
        for lbl in lbls
            vals = [value(ffr.kratios, lbl, missing) for ffr in ffrs]
            insertcols!(res, ncol(res) + 1, Symbol(repr(lbl)) => vals)
            if withUnc
                unc = [σ(ffr.kratios, lbl, missing) for ffr in ffrs]
                insertcols!(res, ncol(res) + 1, Symbol('Δ' * repr(lbl)) => unc)
            end
        end
    end
    return res
end

function NeXLUncertainties.asa(::Type{DataFrame}, fr::FitResult; withUnc = false)
    rois = labels(fr)
    res =  DataFrame(
        ROI = rois,
        k = map(l -> value(fr.kratios, l), rois),
    )
    if withUnc
        res[:, "σ(k)"] = map(l -> σ(fr.kratios, l), rois)
    end
    return res
end

function DataAPI.describe(ffrs::Vector{<:FitResult})::DataFrame
    df = asa(DataFrame, ffrs)[:, 2:end]
    desc = describe(df, :mean, :std, :min, :q25, :median, :q75, :max)
    lbls = filter(
        lbl -> lbl isa CharXRayLabel,
        sort(collect(Set(mapreduce(labels, append!, ffrs)))),
    )
    insertcols!(desc, 4, :hetero => heterogeneity.(Ref(ffrs), lbls))
    return desc
end

"""
    FilterFitResult

Represents the result of fitting a FilteredUnknownW to a FilteredUnknown.

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
    peakback::Dict{ReferenceLabel,NTuple{3,Float64}}
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

"""
    spectrum(ffr::FilterFitResult)::Spectrum

Returns the original unknown spectrum.
"""
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
    peak, back, _ = ffr.peakback[klabel]
    return (peak * (energy(klabel.roi.stop, unk) - energy(klabel.roi.start, unk))) /
           (backwidth * back)
end


"""
    NeXLUncertainties.asa(
        ::Type{DataFrame}, 
        ffr::FilterFitResult; 
        charOnly::Bool=true, 
        material=nothing,
        mc = XPP, fc=ReedFluorescence,
        columns = () # Selected from ( :roi, :peakback, :counts )
    )::DataFrame

Tabulate details about each region-of-interest in the 'FilterFitResult' in a 'DataFrame'.
  * If charOnly then only display characteristic X-ray data (not escapes etc.)
  * If `material` is a Material then the computed k-ratio (KCalc) will also be tabulated along with kmeas/kcalc (KoKCalc).
  * columns - :roi - ROI for peak (Start, Stop)
              :peakback - characteristic intensity, background intensity and PtoB over a 10 eV/channel (Peak, Back, PtoB)
              :counts - Characteristic intensity in peak (Counts)
              :dose - :LiveTime, :ProbeCurrent, :Dose properties from the spectrum
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    ffr::FilterFitResult;
    charOnly::Bool = true,
    material::Union{Material,Nothing} = nothing,
    columns::Tuple = (), # ( :roi, :peakback, :counts, :dose)
    mc = XPP, fc=ReedFluorescence,
)::DataFrame
    sl = NeXLUncertainties.sortedlabels(ffr.kratios)
    unkspec = ffr.label.spectrum
    if charOnly
        sl = filter(lbl->lbl isa CharXRayLabel, sl)
    end
    res = DataFrame(
        Spectrum = [ ffr.label for _ in sl ],
        Feature = sl,
        Reference = [properties(lbl)[:Name] for lbl in sl],
        K = [ NeXLUncertainties.value(ffr.kratios, lbl) for lbl in sl],
        dK = [ NeXLUncertainties.σ(ffr.kratios, lbl) for lbl in sl]
    ) 
    if :roi in columns
        insertcols!(res, 4, :Start => [ lbl.roi.start for lbl in sl])
        insertcols!(res, 5, :Stop => [ lbl.roi.stop for lbl in sl])
    end
    if :peakback in columns
        insertcols!(res, :Peak => [ ffr.peakback[lbl][1]-ffr.peakback[lbl][2] for lbl in sl] )
        insertcols!(res, :Back => [ ffr.peakback[lbl][2] for lbl in sl] )
        insertcols!(res, :PtoB => [ peaktobackground(ffr, lbl) for lbl in sl] )
    end
    if material isa Material
        function kcalc(lbl) 
            stdspec = lbl.spectrum
            zcu = zafcorrection(
                mc,
                fc,
                NullCoating,
                material,
                lbl.xrays,
                unkspec[:BeamEnergy],
            )
            zcs = zafcorrection(
                mc,
                fc,
                NullCoating,
                stdspec[:Composition],
                lbl.xrays,
                stdspec[:BeamEnergy],
            )
            k(zcu, zcs, unkspec[:TakeOffAngle], stdspec[:TakeOffAngle])
        end
        insertcols!(res, :KCalc => kcalc.(sl))
        insertcols!(res, :KoKcalc => [ r[:K]/r[:KCalc] for r in eachrow(res) ])
    end
    if :dose in columns
        dt(spec) = 100.0*(spec[:RealTime]-spec[:LiveTime])/spec[:RealTime]
        insertcols!(res, 2, :LiveTime => [ get(unkspec, :LiveTime, missing) for _ in sl ])
        insertcols!(res, 3, :ProbeCurrent => [ get(unkspec, :ProbeCurrent, missing) for _ in sl ])
        insertcols!(res, 4, :DeadPct => [ dt(unkspec) for lbl in sl ])

    end
    if :counts in columns
        insertcols!(res, :Counts => [ ffr.peakback[lbl][1]-ffr.peakback[lbl][2] for lbl in sl] )
        insertcols!(res, :RefCountsPernAs => [ ffr.peakback[lbl][3]/(lbl.spectrum[:ProbeCurrent]*lbl.spectrum[:LiveTime]) for lbl in sl] )
        insertcols!(res, :CountsPernAs => [ (ffr.peakback[lbl][1]-ffr.peakback[lbl][2]) / (lbl.spectrum[:ProbeCurrent]*lbl.spectrum[:LiveTime]) for lbl in sl ])
    end
    return res
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
    properties(fr.label.spectrum)[:Coating] = estimatecoating(substrate, coating, kratio(fr, cxr), mc)
end