# Precompiling NeXLSpectrum

Unfortunately, when you first start a Julia REPL (or Jupyter/Pluto notebook or Weave document) to work with NeXLSpectrum, 
Julia needs to load and compile many libraries.  This can be tediously slow particularly if you add Gadfly and DataFrames
into the mix.

To mitigate this, it is useful to using PackageCompiler to create a sysimage containing compiled versions of
NeXLUncertainties, NeXLCore, NeXLSpectrum, Gadfly and DataFrames.

This process takes time (about 10 minutes on my laptop) but makes Julia much more responsive.

Here is the secret incantation.
```julia
using PackageCompiler, NeXLSpectrum
PackageCompiler.create_sysimage(
    [ "Gadfly", "DataFrames", "BoteSalvatICX", "FFAST", "NeXLUncertainties", "NeXLCore", "NeXLMatrixCorrection", "NeXLSpectrum" ]; 
    sysimage_path=joinpath(homedir(), ".julia", "NeXLSysimage.dll"),
    precompile_execution_file=joinpath(pkgdir(NeXLSpectrum),"scripts", "precompile.jl"))
```
where `precompile.jl` provides additional usage data to help ensure that the right functions are included.

But this alone isn't enough.  Now you need to tell Julia that you want to load this sysimage.  For this you need to add
the `-J` argument when you invoke julia.  Something like this
```
> julia -J/$HOME$/.julia/NeXLSysmage.dll
```
where `#HOME` is the path in which `.julia` is located.

It is worthwhile to understand the [limitations](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html#Drawbacks-to-custom-sysimages) associated with PackageCompiler and sysimages.