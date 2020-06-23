using Documenter, NeXLSpectrum

let curr=pwd()
    try
        include(joinpath("..","weave","buildweave.jl"))
    finally
        cd(curr)
    end
end

makedocs(
    modules = [NeXLSpectrum],
    sitename = "NeXLSpectrum.jl",
    pages = [
        "Home" => "index.md",
        "Fitting K412" => "K412fit.md",
        "Fitting K412 (quick fit)" => "K412quick.md",
        #"Quantifying AMD-6005a glass" => "quantAMDglass.md",
        "Fitting XRF Spectra" => "XRFSpectra.md",
        "Lovely Error Bars" => "errorbars.md",
        "Modeling the Continuum" => "continuummodel.md"
     ]
)

map(name->rm(joinpath("src","$(splitext(name)[1]).md")), names)

# deploydocs(repo = "github.com/NeXLSpectrum.jl.git")
