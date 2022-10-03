"""
    NeXLMatrixCorrection.quantify(
        ffr::FitResult,
        iteration::Iteration = Iteration(mc=XPP, fc=ReedFluorescence, cc=Coating);
        strip::AbstractVector{Element} = Element[],
        kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
        coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
    )::IterationResult

Facilitates converting `FilterFitResult` or `BasicFitResult` objects into estimates of composition by extracting
k-ratios from measured spectra and applying matrix correction algorithms.
"""
function NeXLMatrixCorrection.quantify(
    ffr::FitResult;
    strip::AbstractVector{Element} = Element[],
    iteration::Iteration = Iteration(),
    kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
    coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
)::IterationResult
    krs = optimizeks(kro, filter(kr -> !(element(kr) in strip), kratios(ffr)))
    return quantify(ffr.label, krs, iteration, coating=coating)
end

"""
    NeXLMatrixCorrection.quantify(
        spec::Union{Spectrum,AbstractVector{<:Spectrum}},
        ffp::FilterFitPacket;
        strip::AbstractVector{Element} = Element[],
        iteration::Iteration = Iteration(mc=XPP, fc=ReedFluorescence, cc=Coating),
        kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
        coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
    )::IterationResult

Failitates quantifying spectra.  First filter fits and then matrix corrects.
"""
NeXLMatrixCorrection.quantify(
    spec::Spectrum,
    ffp::FilterFitPacket;
    kwargs...,
) = quantify(fit_spectrum(spec, ffp); kwargs...)

NeXLMatrixCorrection.quantify(
    specs::AbstractVector{<:Spectrum},
    ffp::FilterFitPacket;
    kwargs...,
)::Vector{IterationResult} = #
    map(spec->quantify(fit_spectrum(spec, ffp); kwargs...), specs)