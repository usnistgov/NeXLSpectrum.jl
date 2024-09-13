module NeXLSpectrumDataFrameExt

import NeXLUncertainties
import NeXLCore
import NeXLMatrixCorrection
import NeXLSpectrum
import DataFrames
import PeriodicTable
import DataAPI


"""
Special `DataFrames.DataFrame` support for `NeXLSpectrum` types.

    DataFrames.DataFrame(spec::Spectrum)

    DataFrames.DataFrame(::Type{Spectrum}, spec::AbstractArray{Spectrum})::DataFrame

    DataFrames.DataFrame(spec::Spectrum; properties::Bool = false)
    
    DataFrames.DataFrame(::Type{Element}, spec::AbstractDict{Element, Spectrum})::DataFrame

    DataFrames.DataFrame(ffp::FilterFitPacket)
    
    DataFrames.DataFrame(fr::FitResult; withUnc = false)
    
    DataFrames.DataFrame(ffp::DirectReferences)
"""
function DataFrames.DataFrame(spec::NeXLSpectrum.Spectrum{T}) where {T}
    df = DataFrames.DataFrame(:Property => String[], :Value => String[])
    push!(df, ("Name", "$(spec[:Name])"))
    push!(df, ("Calibration", "$(spec.energy)"))
    push!(df, ("Beam energy", "$(get(spec, :BeamEnergy, missing)/1000.0) keV"))
    push!(df, ("Probe current", "$(get(spec, :ProbeCurrent, missing)) nA"))
    push!(df, ("Live time", "$(get(spec, :LiveTime, missing)) s"))
    push!(df, ("Coating", "$(get(spec,:Coating, nothing))"))
    push!(df, ("Detector", "$(get(spec, :Detector, missing))"))
    push!(df, ("Comment", "$(get(spec, :Comment, missing))"))
    push!(df, ("Integral", "$(NeXLSpectrum.integrate(spec)) counts"))
    comp = get(spec, :Composition, missing)
    if !ismissing(comp)
        push!(df, ("Composition", "$(comp)"))
        det = get(spec, :Detector, missing)
        if !ismissing(det)
            coating = get(spec, :Coating, missing)
            comp2 = Set(keys(comp))
            if !ismissing(coating)
                union!(comp2, keys(coating.material))
            end
            for elm1 in comp2
                for ext1 in NeXLSpectrum.extents(elm1, det, 1.0e-4)
                    intersects = []
                    for elm2 in comp2
                        if elm2 ≠ elm1
                            for ext2 in NeXLSpectrum.extents(elm2, det, 1.0e-4)
                                if length(intersect(ext1, ext2)) > 0
                                    push!(intersects, "$(elm2.symbol)[$(ext2)]")
                                end
                            end
                        end
                    end
                    if length(intersects) > 0
                        push!(df, ("ROI $(elm1.symbol)[$(ext1)]","Intersects $(join(intersects,", "))"))
                    else
                        p, b = NeXLSpectrum.peak(spec, ext1), NeXLSpectrum.background(spec, ext1)
                        σ = p / sqrt(b)
                        push!(df, ("ROI $(elm1.symbol)[$(ext1)]", "$(round(Int,p)) counts over $(round(Int,b)) counts - σ = $(round(Int,σ))"))
                    end
                end
            end
        end
    end
    return df
end

"""
    DataFrames.DataFrame(
        ffr::FilterFitResult; 
        charOnly::Bool=true, 
        material=nothing,
        mc = XPP, fc=ReedFluorescence,
        columns = ( :counts, ) # Selected from ( :roi, :peakback, :counts, :dose )
    )::DataFrames.DataFrame

Tabulate details about each region-of-interest in the 'FilterFitResult' in a 'DataFrame'.
  * If charOnly then only display characteristic X-ray data (not escapes etc.)
  * If `material` is a Material then the computed k-ratio (KCalc) will also be tabulated along with kmeas/kcalc (KoKCalc).
  * columns - Select optional column outputs (see below)

Columns:

  * Spectrum : UnknownLabel - Identifies the fit spectrum
  * Feature  : Label - Identifies the fit feature (Vector{CharXRay} or other)
  * Reference: String - Name of reference against which :Spectrum was fit over :Feature
  * K : Float64 - The multiplicative fit constant
  * dK : Float64 - The 1σ uncertainty in :K
  * :Start : Int - Start index for fit channels (:roi ∈ columns)
  * :Stop : Int - Stop index for fit channels  (:roi ∈ columns)
  * :Counts : Float64 - Total counts in characteristic peak (peak-back) (:peakback || :counts ∈ columns)
  * :Back : Float64 - Total counts in background under the characteristic peak (:peakback ∈ columns)
  * :PtoB : Float64 - Peak-to-Background assuming 10 eV/channel (:peakback ∈ columns)
  * :KCalc : Float64 - Computed k-ratio assuming a composition. (Requires `material` argument to be specified.)
  * :KoKcalc : Float64 - Ratio of measured/computed k-ratio.  (Requires `material` argument to be specified.)
  * :LiveTime : Float64 - Acquisiton live time (s) (:dose ∈ columns)
  * :ProbeCurrent : Float64 - Probe current (nA) (:dose ∈ columns)
  * :DeadPct : Float64 - Dead time in ProbeCurrent (:dose ∈ columns)
  * :RefCountsPernAs : Float64 - Estimated counts in :Reference in :Feature per unit dose.  (:counts ∈ columns)
  * :CountsPernAs : Float64 - Estimated counts in :Spectrum in :Feature per unit dose.  (:counts ∈ columns)

  Note: PtoB is defined as "The peak-to-background ratio as determined from the raw and residual spectra integrated over
  the fit region-of-interest and scaled to 10 eV of continuum." Why? Cause.
"""
function DataFrames.DataFrame(
    ffr::NeXLSpectrum.FilterFitResult;
    charOnly::Bool = true,
    material::Union{NeXLCore.Material,Nothing} = nothing,
    columns::Tuple = ( :counts, ), # ( :roi, :peakback, :counts, :dose)
    mc = NeXLMatrixCorrection.XPP, fc=NeXLMatrixCorrection.ReedFluorescence,
)::DataFrames.DataFrame
    sl = charOnly ? #
        filter(lbl->lbl isa NeXLSpectrum.CharXRayLabel, NeXLUncertainties.sortedlabels(ffr.kratios)) : #
        NeXLUncertainties.sortedlabels(ffr.kratios)
    unkspec = ffr.label.spectrum
    # Build the base result DataFrame
    res = DataFrames.DataFrame(
        Spectrum = [ ffr.label for _ in sl ],
        Feature = sl,
        Reference = [NeXLCore.properties(lbl)[:Name] for lbl in sl],
        K = [ NeXLUncertainties.value(ffr.kratios, lbl) for lbl in sl],
        dK = [ NeXLUncertainties.σ(ffr.kratios, lbl) for lbl in sl]
    ) 
    # Add optional columns
    if :roi ∈ columns
        DataFrames.insertcols!(res, 4, :Start => [ lbl.roi.start for lbl in sl])
        DataFrames.insertcols!(res, 5, :Stop => [ lbl.roi.stop for lbl in sl])
    end
    if :peakback ∈ columns
        pb = ffr.peakback()
        DataFrames.insertcols!(res, :Counts => [ (pb[lbl][1]-pb[lbl][2]) for lbl in sl] )
        DataFrames.insertcols!(res, :Back => [ pb[lbl][2] for lbl in sl] )
        DataFrames.insertcols!(res, :PtoB => [ NeXLSpectrum.peaktobackground(ffr, lbl) for lbl in sl] )
    end
    if material isa NeXLCore.Material
        function kcalc(lbl) 
            stdspec = lbl.spectrum
            zc = NeXLMatrixCorrection.zafcorrection(
                mc,
                fc,
                NeXLMatrixCorrection.NullCoating,
                material,
                stdspec[:Composition],
                lbl.xrays,
                unkspec[:BeamEnergy],
            )
            NeXLMatrixCorrection.k(zc..., unkspec[:TakeOffAngle], stdspec[:TakeOffAngle])
        end
        DataFrames.insertcols!(res, :KCalc => kcalc.(sl))
        DataFrames.insertcols!(res, :KoKcalc => [ r[:K]/r[:KCalc] for r in eachrow(res) ])
    end
    if :dose ∈ columns
        DataFrames.insertcols!(res, 2, :LiveTime => [ get(unkspec, :LiveTime, missing) for _ in sl ])
        DataFrames.insertcols!(res, 3, :RealTime => [ get(unkspec, :RealTime, missing) for _ in sl ])
        DataFrames.insertcols!(res, 4, :ProbeCurrent => [ get(unkspec, :ProbeCurrent, missing) for _ in sl ])
    end
    if :counts ∈ columns
        pb = ffr.peakback()
        if !(:peakback ∈ columns)
            DataFrames.insertcols!(res, :Counts => [ pb[lbl][1]-pb[lbl][2] for lbl in sl] )
        end
        DataFrames.insertcols!(res, :RefCountsPernAs => [ pb[lbl][3] for lbl in sl] )
        DataFrames.insertcols!(res, :CountsPernAs => [ (pb[lbl][1]-pb[lbl][2]) / NeXLSpectrum.dose(unkspec) for lbl in sl ])
    end
    return res
end

function DataFrames.DataFrame(ffp::NeXLSpectrum.FilterFitPacket)
    # charonly over roi, data over ffroi
    croi(fref) = max(1,first(fref.roi)-first(fref.ffroi)+1):min(length(fref.data), last(fref.roi)-first(fref.ffroi)+1)
    function p2b(fref)
        @assert length(croi(fref))==length(fref.roi)
        sum(fref.charonly) / (sum(fref.data[croi(fref)])-sum(fref.charonly))
    end 
    function s2n(fref)
        @assert length(croi(fref))==length(fref.roi)
        sum(fref.charonly) / sqrt(sum(fref.data[croi(fref)])-sum(fref.charonly))
    end 
    DataFrames.DataFrame(
        Symbol("Spectrum") => [ NeXLSpectrum.name(fr.label.spectrum) for fr in ffp.references ],
        Symbol("Beam Energy (keV)") => [ get(fr.label.spectrum, :BeamEnergy, missing)/1000.0 for fr in ffp.references ],
        Symbol("Probe Current (nA)") => [ get(fr.label.spectrum, :ProbeCurrent, missing) for fr in ffp.references ],
        Symbol("Live Time (s)") => [ get(fr.label.spectrum, :LiveTime, missing) for fr in ffp.references ],
        Symbol("Material") => [ get(fr.label.spectrum, :Composition, nothing) for fr in ffp.references],
        Symbol("Lines") => [ fr.label.xrays for fr in ffp.references],
        Symbol("Coating") => [ get(fr.label.spectrum, :Coating, nothing) for fr in ffp.references ],
        Symbol("ROI") => [ fr.roi for fr in ffp.references],
        Symbol("Full ROI") => [ fr.ffroi for fr in ffp.references],
        Symbol("P-to-B") => p2b.(ffp.references), 
        Symbol("S-to-N") => s2n.(ffp.references) 
    )
end

function DataFrames.DataFrame(ffp::NeXLSpectrum.DirectReferences)
    DataFrames.DataFrame(
        :Type => [ "$(fr.label)" for fr in ffp.references ],
        :Spectrum => [ NeXLSpectrum.name(fr.label.spectrum) for fr in ffp.references ],
        Symbol("Beam Energy (keV)") => [ get(fr.label.spectrum, :BeamEnergy, missing)/1000.0 for fr in ffp.references ],
        Symbol("Probe Current (nA)") => [ get(fr.label.spectrum, :ProbeCurrent, missing) for fr in ffp.references ],
        Symbol("Live Time (s)") => [ get(fr.label.spectrum, :LiveTime, missing) for fr in ffp.references ],
        :Material => [ get(fr.label.spectrum, :Composition, nothing) for fr in ffp.references],
        :Lines => [ fr.label.xrays for fr in ffp.references],
        :ROI => [ fr.roi for fr in ffp.references],
    )
end

function DataAPI.describe(ffrs::Vector{Union{<:NeXLSpectrum.FitResult,<:NeXLSpectrum.FilterFitResult}})::DataFrames.DataFrame
    df = DataFrames.DataFrame(ffrs)[:, 2:end]
    desc = DataAPI.describe(df, :mean, :std, :min, :q25, :median, :q75, :max)
    lbls = filter(
        lbl -> lbl isa NeXLSpectrum.CharXRayLabel,
        sort(collect(Set(mapreduce(NeXLSpectrum.labels, append!, ffrs)))),
    )
    DataFrames.insertcols!(desc, 4, :hetero => heterogeneity.(Ref(ffrs), lbls))
    return desc
end

function DataFrames.DataFrame(fr::NeXLSpectrum.FitResult; withUnc = false)
    rois = NeXLSpectrum.labels(fr)
    res =  DataFrames.DataFrame(
        Spectrum = fill(fr.label, length(rois)),
        Feature = rois,
        K = map(l -> NeXLUncertainties.value(fr.kratios, l), rois),
    )
    withUnc && DataFrames.insertcols!(res, :dK => map(l -> NeXLUncertainties.σ(fr.kratios, l), rois) )
    return res
end


"""
    DataFrames.DataFrame(
        ::Type{NeXLSpectrum.FitResult},
        ffrs::AbstractVector{<:NeXLSpectrum.FitResult};
        charOnly = true,
        withUnc = false,
        format = :normal # :pivot or :long
    )::DataFrames.DataFrame

Tabulate a generic `FilterFitResult`, `DirectFitResult` or other `FitResult` as a DataFrame.

Format:
  * :normal - One row per spectrum, one column per k-ratio
  * :pivot - One row per ROI, one column per spectrum (optional: column for 1σ uncertainty on k-ratio)
  * :long - One row per spectrum, feature, measured k-ratio

"""
function DataFrames.DataFrame(
    ::Type{NeXLSpectrum.FitResult},
    ffrs::AbstractVector{<:NeXLSpectrum.FitResult};
    charOnly = true,
    withUnc = false,
    format = :normal # :pivot or :long
)::DataFrames.DataFrame
    lbls = sort(unique(mapreduce(NeXLUncertainties.labels, append!, ffrs)))
    lbls = charOnly ? filter(lbl -> lbl isa NeXLSpectrum.CharXRayLabel, lbls) : lbls
    if format==:pivot
        res = DataFrames.DataFrame(ROI = lbls)
        for ffr in ffrs
            vals = [ NeXLUncertainties.value(ffr.kratios, lbl, missing) for lbl in lbls ]
            DataFrames.insertcols!(res, DataAPI.ncol(res) + 1, Symbol(repr(ffr.label)) => vals)
            if withUnc
                unc = [NeXLUncertainties.σ(ffr.kratios, lbl, missing) for lbl in lbls]
                DataFrames.insertcols!(res, DataAPI.ncol(res) + 1, Symbol('Δ' * repr(ffr.label)) => unc)
            end
        end
    elseif format==:long    
        if withUnc
            res = DataFrames.DataFrame(Spectrum=String[], ROI = String[], k = Float64[], dk = Float64[]) 
            for ffr in ffrs, lbl in lbls
                push!(res, ( repr(ffr.label), repr(lbl), NeXLUncertainties.value(ffr.kratios, lbl, missing), NeXLUncertainties.σ(ffr.kratios, lbl, missing) ))
            end
        else
            res = DataFrames.DataFrame(Spectrum=String[], ROI = String[], k = Float64[]) 
            for ffr in ffrs, lbl in lbls
                push!(res, ( repr(ffr.label), repr(lbl), NeXLUncertainties.value(ffr.kratios, lbl, missing)))
            end
        end
    else
        rowLbls = [repr(ffr.label) for ffr in ffrs]
        res = DataFrames.DataFrame(Symbol("Spectra") => rowLbls)
        for lbl in lbls
            vals = [NeXLUncertainties.value(ffr.kratios, lbl, missing) for ffr in ffrs]
            DataFrames.insertcols!(res, DataAPI.ncol(res) + 1, Symbol(repr(lbl)) => vals)
            if withUnc
                unc = [NeXLUncertainties.σ(ffr.kratios, lbl, missing) for ffr in ffrs]
                DataFrames.insertcols!(res, DataAPI.ncol(res) + 1, Symbol('Δ' * repr(lbl)) => unc)
            end
        end
    end
    return res
end
DataFrames.DataFrame(
    ::Type{NeXLSpectrum.FilterFitResult},
    ffrs::AbstractVector{<:NeXLSpectrum.FilterFitResult};
    kwargs...
) = DataFrames.DataFrame(NeXLSpectrum.FitResult, ffrs; kwargs...)

function DataFrames.DataFrame(::Type{NeXLSpectrum.Spectrum}, specs::AbstractVector{<:NeXLSpectrum.Spectrum})::DataFrames.DataFrame
    _asname(comp) = ismissing(comp) ? missing : NeXLCore.name(comp)
    unf, unl, uns = Union{Float64,Missing}, Union{NeXLCore.Film,Nothing}, Union{String,Missing}
    nme, e0, pc, lt, rt, coat, integ, comp =
        String[], unf[], unf[], unf[], unf[], unl[], Float64[], uns[]
    for spec in specs
        push!(nme, spec[:Name])
        push!(e0, get(spec, :BeamEnergy, missing))
        push!(pc, get(spec, :ProbeCurrent, missing))
        push!(lt, get(spec, :LiveTime, missing))
        push!(rt, get(spec, :RealTime, missing))
        push!(coat, get(spec, :Coating, nothing))
        push!(integ, NeXLSpectrum.integrate(spec))
        push!(comp, _asname(get(spec, :Composition, missing)))
    end
    return DataFrames.DataFrame(
        Name = nme,
        BeamEnergy = e0,
        ProbeCurrent = pc,
        LiveTime = lt,
        RealTime = rt,
        Coating = coat,
        Integral = integ,
        Material = comp,
    )
end


function DataFrames.DataFrame(::Type{NeXLSpectrum.Spectrum}, stds::AbstractDict{PeriodicTable.Element,<:NeXLSpectrum.Spectrum})::DataFrames.DataFrame
    _asname(comp) = ismissing(comp) ? missing : NeXLCore.name(comp)
    elms = sort(collect(keys(stds)))
    return DataFrames.DataFrame(
        Element = map(el->el.symbol, elms),
        Z = map(el->el.number, elms),
        Name = map(el->stds[el][:Name], elms),
        Material = map(el->_asname(get(stds[el], :Composition, missing)), elms),
        MassFrac = map(el->ismissing(stds[el][:Composition]) ? missing : stds[el][:Composition][el], elms),
        BeamEnergy = map(el->get(stds[el], :BeamEnergy, missing), elms),
        ProbeCurrent = map(el->get(stds[el], :ProbeCurrent, missing), elms),
        LiveTime = map(el->get(stds[el], :LiveTime, missing), elms),
        RealTime = map(el->get(stds[el], :RealTime, missing), elms),
        Coating = map(el->get(stds[el], :Coating, nothing), elms),
        Integral = map(el->NeXLSpectrum.integrate(stds[el]), elms),
    )
end

function NeXLSpectrum.suitability(elm::PeriodicTable.Element, det::NeXLSpectrum.EDSDetector; maxE=30.0e3, minC=0.1, latex=false)
    mats = NeXLMatrixCorrection.getstandards(elm, minC)
    wrap(matname, ltx) = ltx ? "\\ce{$matname}" : matname
    mark(b) = b ? (latex ? "\\checkmark" : "✓") : (latex ? "\\xmark" : "✗")
    mats = filter(m->NeXLUncertainties.value(m[elm])>0.0, mats)
    # Find all elemental ROIs
    rois = Dict{Vector{NeXLCore.CharXRay}, Vector{NeXLCore.Material}}()
    for (cxrs, _) in NeXLSpectrum.suitablefor(elm, Set( (elm, ) ), det, maxE=maxE, ampl=1.0e-5, warnme=false)
      rois[cxrs] = NeXLCore.Material[]
    end
    for mat in mats
      for (cxrs, _) in NeXLSpectrum.suitablefor(elm, mat, det, maxE=maxE, ampl=1.0e-5, warnme=false)
        push!(rois[cxrs], mat) 
      end
    end
    cxrss = collect(keys(rois))
    sort!(cxrss, lt = (a,b) -> NeXLCore.energy(NeXLCore.brightest(a))<NeXLCore.energy(NeXLCore.brightest(b)))
    lmats = collect(mats)
    res = DataFrames.DataFrame(
        :Material=>map(m->wrap(NeXLCore.name(m), latex), lmats),
        :MassFrac =>map(m->round(NeXLUncertainties.value(m[elm]), digits=3), lmats), 
        :Count=>map(m->count(cxrs->m in rois[cxrs], cxrss), lmats), 
        map(cxrs -> Symbol(cxrs)=>map(m-> mark(m in rois[cxrs]), lmats), cxrss)...
    )
    sort!(res, [ :Count, :MassFrac ], rev=true)
end
function NeXLSpectrum.suitability(elm::PeriodicTable.Element, mats::Union{AbstractSet{<:NeXLCore.Material}, AbstractVector{<:NeXLCore.Material}}, det::NeXLSpectrum.EDSDetector; maxE=30.0e3, latex=false)
    wrap(matname, ltx) = ltx ? "\\ce{$matname}" : matname
    mark(b) = b ? (latex ? "\\checkmark" : "✓") : (latex ? "\\xmark" : "✗")
    mats = filter(m->NeXLUncertainties.value(m[elm])>0.0, mats)
    # Find all elemental ROIs
    rois = Dict{Vector{NeXLCore.CharXRay}, Vector{NeXLCore.Material}}()
    for (cxrs, _) in NeXLSpectrum.suitablefor(elm, Set( (elm, ) ), det, maxE=maxE, ampl=1.0e-5, warnme=false)
      rois[cxrs] = NeXLCore.Material[]
    end
    for mat in mats
      for (cxrs, _) in NeXLSpectrum.suitablefor(elm, mat, det, maxE=maxE, ampl=1.0e-5, warnme=false)
        push!(rois[cxrs], mat) 
      end
    end
    cxrss = collect(keys(rois))
    sort!(cxrss, lt = (a,b) -> NeXLCore.energy(NeXLCore.brightest(a))<NeXLCore.energy(NeXLCore.brightest(b)))
    lmats = collect(mats)
    res = DataFrames.DataFrame(
        :Material=>map(m->wrap(NeXLCore.name(m), latex), lmats),
        :MassFrac =>map(m->round(NeXLUncertainties.value(m[elm]), digits=3), lmats), 
        :Count=>map(m->count(cxrs->m in rois[cxrs], cxrss), lmats), 
        map(cxrs -> Symbol(cxrs)=>map(m-> mark(m in rois[cxrs]), lmats), cxrss)...
    )
    sort!(res, [ :Count, :MassFrac ], rev=true)
end

NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, ffp::NeXLSpectrum.DirectReferences) = DataFrames.DataFrame(ffp)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, ffrs::AbstractVector{<:NeXLSpectrum.FitResult}; kwargs...) = DataFrames.DataFrame(NeXLSpectrum.FitResult, ffrs; kwargs...)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, fr::NeXLSpectrum.FitResult; kwargs...) = DataFrames.DataFrame(fr; kwargs...)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, ffp::NeXLSpectrum.FilterFitPacket) = DataFrames.DataFrame(ffp)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, specs::AbstractArray{<:NeXLSpectrum.Spectrum}) = DataFrames.DataFrame(NeXLSpectrum.Spectrum, specs)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, spec::NeXLSpectrum.Spectrum; kwargs...) = DataFrames.DataFrame(spec; kwargs...)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, ffr::NeXLSpectrum.FilterFitResult; kwargs...) = DataFrames.DataFrame(ffr; kwargs...)
NeXLUncertainties.asa(::Type{DataFrames.DataFrame}, stds::AbstractDict{PeriodicTable.Element,NeXLSpectrum.Spectrum}) = DataFrames.DataFrame(NeXLSpectrum.Spectrum, stds)    

end