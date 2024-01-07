module NeXLSpectrumUnicodePlotsExt

using NeXLSpectrum
using UnicodePlots

function UnicodePlots.lineplot(sps::Spectrum...; kvargs...)
    emax = maximum(s->get(s, :BeamEnergy, 20.0e3), sps)
    xlim = :xlim in keys(kvargs) ? kvargs[:xlim] : (0.0, emax)
    plt = nothing
    for sp in sps
        plt = isnothing(plt) ?
            lineplot(energyscale(sp), counts(sp); kvargs..., xlim=xlim) :
            lineplot!(plt, energyscale(sp), counts(sp))
    end
    return plt
end

end