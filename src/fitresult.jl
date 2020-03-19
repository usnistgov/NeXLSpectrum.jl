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
struct FilterFitResult
    label::UnknownLabel
    kratios::UncertainValues
    roi::UnitRange{Int}
    raw::Vector{Float64}
    residual::Vector{Float64}
    peakback::Dict{<:ReferenceLabel,NTuple{2,Float64}}
end

"""
    kratios(ffr::FilterFitResult)

The k-ratios as a UncertainValues object
"""
kratios(ffr::FilterFitResult)::UncertainValues = ffr.kratios
unknown(ffr::FilterFitResult)::UnknownLabel = ffr.label
"""
    residual(ffr::FilterFitResult)::Spectrum

A Spectrum containing the histogram representing the unknown spectrum
minus the fitted characteristic peaks shapes times the best fit coefficient.
"""
function residual(ffr::FilterFitResult)::Spectrum
    props = copy(ffr.label.spec.properties)
    props[:Name]="Residual[$(ffr.label.spec.properties[:Name])]"
    return Spectrum(ffr.label.spec.energy, ffr.residual, props)
end

"""
    characteristiccounts(ffr::FiterFitResult, strip)

Number of spectrum counts that were accounted for by the fitted elements with the `strip` Element(s) removed.
"""
characteristiccounts(ffr::FilterFitResult,strip) = sum(element(ref) in strip ? 0.0 : v[1]-v[2] for (ref,v) in ffr.peakback) # sum(ffr.raw[ffr.roi]-ffr.residual[ffr.roi])

"""
    peaktobackground(ffr::FilterFitResult, backwidth::Float64=10.0)::Float64

The peak-to-background ratio as determined from the raw and residual spectra integrated over
the fit region-of-interest and scaled to <code>backwidth</code> eV of continuum (nominally 10 eV).
"""
function peaktobackground(ffr::FilterFitResult, klabel::ReferenceLabel, backwidth::Float64 = 10.0)::Float64
    unk = spectrum(unknown(ffr))
    peak, back = ffr.peakback[klabel]
    return (peak * (energy(klabel.roi.stop, unk) - energy(klabel.roi.start, unk))) / (backwidth * back)
end

NeXLUncertainties.asa(::Type{DataFrame}, ffrs::Array{FilterFitResult}, withUnc = false)::DataFrame =
    hcat(DataFrame(Name = [r.label for r in ffrs]), asa(DataFrame, [r.kratios for r in ffrs], withUnc))

function NeXLUncertainties.asa(::Type{DataFrame}, ffr::FilterFitResult)::DataFrame
    lbl, klbl, std, kr, dkr, roi1, roi2, peak, back, ptob = UnknownLabel[],
        ReferenceLabel[],
        String[],
        Float64[],
        Float64[],
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[]
    for kl in sortedlabels(ffr.kratios)
        push!(lbl, ffr.label)
        push!(std, spectrum(kl)[:Name])
        push!(klbl, kl)
        push!(roi1, kl.roi.start)
        push!(roi2, kl.roi.stop)
        push!(kr, value(kl, ffr.kratios))
        push!(dkr, σ(kl, ffr.kratios))
        pb = ffr.peakback[kl]
        push!(peak, pb[1])
        push!(back, pb[2])
        push!(ptob, peaktobackground(ffr, kl))
    end
    return DataFrame(
        Label = lbl,
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

Base.show(io::IO, ffr::FilterFitResult) = print(io, "$(ffr.label)")

NeXLUncertainties.value(label::ReferenceLabel, ffr::FilterFitResult) = value(label, ffr.kratios)
NeXLUncertainties.σ(label::ReferenceLabel, ffr::FilterFitResult) = σ(label, ffr.kratios)
NeXLUncertainties.getindex(ffr::FilterFitResult, label::ReferenceLabel) = getindex(ffr.kratios, label)
Base.keys(ffr::FilterFitResult) = keys(ffr.kratios)
NeXLUncertainties.labels(ffr::FilterFitResult) = labels(ffr.kratios)

"""
    filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
function filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}
    scaled(ff) = (value(ff.identifier, fit) * (ff.scale / unk.scale)) * extract(ff, unk.ffroi)
    return unk.filtered - mapreduce(scaled, +, ffs)
end
