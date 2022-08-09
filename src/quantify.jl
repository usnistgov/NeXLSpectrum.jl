"""
    NeXLMatrixCorrection.quantify(
        ffr::FitResult;
        strip::AbstractVector{Element} = Element[],
        iteration::Iteration = Iteration(mc=XPP, fc=ReedFluorescence, cc=Coating),
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
function NeXLMatrixCorrection.quantify(
    spec::Spectrum,
    ffp::FilterFitPacket,
    iteration::Iteration = Iteration();
    kwargs...,
)::IterationResult
    return quantify(fit_spectrum(spec, ffp); kwargs...)
end
NeXLMatrixCorrection.quantify(
    specs::AbstractVector{<:Spectrum},
    ffp::FilterFitPacket,
    iteration::Iteration = Iteration();
    kwargs...,
)::Vector{IterationResult} = #
    map(spec->quantify(fit_spectrum(spec, ffp), iteration; kwargs...), specs)