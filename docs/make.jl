using Documenter, NeXLSpectrum

makedocs(modules = [NeXLSpectrum], sitename = "NeXLSpectrum.jl")
include(joinpath(@__DIR__,"..","weave","buildweave.jl"))

# deploydocs(repo = "github.com/NeXLSpectrum.jl.git")
