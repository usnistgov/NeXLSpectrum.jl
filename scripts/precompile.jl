using DataFrames
using Gadfly
using BoteSalvatICX
using FFAST
using NeXLUncertainties
using NeXLCore
using NeXLMatrixCorrection
using NeXLSpectrum

include(joinpath(pkgdir(NeXLSpectrum), "test", "runtests.jl"))

specs = [ loadspectrum(joinpath(pkgdir(NeXLSpectrum), "test", "K412 spectra", "III-E K412[$i][4].msa")) for i in 0:4 ]
plot(specs..., autoklms=true, klms = [ n"C", ], edges = [ n"Fe K", n"Ca K" ], coincidences = [ n"O K-L3", ]) |> SVG(tempname(), 10inch, 3inch)

#=
using PackageCompiler, NeXLSpectrum
PackageCompiler.create_sysimage(
    [ "Gadfly", "DataFrames", "BoteSalvatICX", "FFAST", "NeXLUncertainties", "NeXLCore", "NeXLMatrixCorrection", "NeXLSpectrum" ]; 
    sysimage_path=joinpath(homedir(), ".julia", "NeXLSysimage.dll"),
    precompile_execution_file=joinpath(pkgdir(NeXLSpectrum),"scripts", "precompile.jl"))
=#

## To use the sysimage, start Julia with the `-J` option
# julia -JC:\Users\username\.julia\NeXLSysimage.dll"