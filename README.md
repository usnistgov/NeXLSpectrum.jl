# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis
| **Documentation**                        | **Build Status**                  |
|:----------------------------------------:|:---------------------------------:|
| [![][docs-stable-img]][docs-stable-url]  | [![][travis-img]][travis-url]     |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://pages.nist.gov/NeXLSpectrum.jl
[travis-img]: https://travis-ci.com/usnistgov/NeXLSpectrum.jl.svg?branch=master
[travis-url]: https://travis-ci.com/usnistgov/NeXLSpectrum.jl

#### Installation
Install NeXLSpectrum using the Julia package manager
```julia
julia> ]add NeXLSpectrum
```
or

```julia
julia> using Pkg
julia> Pkg.add("NeXLSpectrum")
```

#### Notes
`NeXLSpectrum` is a library of tools for manipulating EDS spectrum within the
NeXL framework. `NeXLSpectrum` depends on `NeXLUncertainties`, `NeXLCore` and
`NeXLMatrixCorrection` and loading `NeXLSpectrum` will also make these
libraries available.

In addition, `NeXLSpectrum` makes extensive use of the third-party libraries
`DataFrames` for data tables and `Gadfly` for plotting.  You will also want to
makes these libraries available.
```julia
julia> ]add DataFrames, Gadfly
```

Primarily, `NeXLSpectrum`
  * Implements the `Spectrum` type to represent individual EDS spectra
    * Reads `Spectrum` objects from disk files (or other streams) in EMSA, Bruker and ASPEX formats
    * Writes `Spectrum` objects to a disk file in EMSA format
  * Provides utilities and other low level tools to interogate and manipulate `Spectrum` objects
  * Implements the `HyperSpectrum` type to represent hyper-spectra (linescan, image, cube, ...)
    * The individual pixels in a hyper-spectrum are visible as `Spectrum` objects
    * Reads `HyperSpectrum` objects from LISPIX-style RPL/RAW files
    * Writes `HyperSpectrum` objects to RPL/RAW files

  * Provides data types to define detector properties
  * Extends `Gadfly.jl` to plot spectra and spectrum-related items
  * Provides algorithms to perform Schamber-style filter-fitting of spectra
    * Implements a basic weighted LLSQ fit algorithm
      * Fits characteristic, escape, Compton and other features
    * Implements a 'vector-based' quick-quant algorithm for processing hyper-spectra

## NeXL and PackageCompiler
Because Julia uses a just-in-time compiler that recompiles libraries each time Julia
is restarted and the `NeXL` libraries are quite large, the first time you perform an 
operation it can be quite slow.  Subsequent calls are much faster.  Fortunately, there 
is a way around this.  The `PackageCompiler` library can be used to build a version of 
Julia in which `NeXLSpectrum` and the libraries on which it depends are build into the 
language much like the base libraries.  I recommend that you use the `PackageCompiler` 
to rebuild Julia each time you update the packages (using `julia> ]up`) or upgrade
to a new version of Julia.

In Windows, this operation looks like:
```julia
julia> using NeXLSpectrum, DataFrames, Gadfly # Load the basic libraries
julia> using PackageCompiler 
julia> PackageCompiler.create_sysimage(
    [ "Gadfly", "DataFrames", "BoteSalvatICX", "FFAST", "NeXLUncertainties", "NeXLCore", "NeXLMatrixCorrection", "NeXLSpectrum" ]; 
    sysimage_path=joinpath(homedir(), ".julia", "NeXLSysimage.dll"),
    precompile_execution_file=joinpath(pkgdir(NeXLSpectrum),"scripts", "precompile.jl"))
```
This operation may take tens-of-minutes to complete but will be worth the effort in terms of saved time
later on.

In other operating systems, you will want to change the extension on the filename in the `sysimage_path`.

This command creates a `sysimage` containing the `NeXL` libraries along with `DataFrames` and `Gadfly`.  You will now
have to tell Julia to use this `sysimage`.   I usually create multiple links to Julia that both use the default `sysimage` 
and the `NeXLSysimage`.  The `sysimage` is selected using the `-J` command line option. In Windows, this looks like:

```
## To use the sysimage, start Julia with the `-J` option
# > julia -J"C:\Users\username\.julia\NeXLSysimage.dll"
```

You can edit the properties of the link used to start Julia to include the `-J` option.  While you are at it,
you might also want to add the `-t` option to start Julia with multiple threads.  For example,
```
## To use the sysimage, start Julia with the `-J` option
# > julia -J"C:\Users\username\.julia\NeXLSysimage.dll" -t 4
```
will start Julia with the `sysimage` and 4 threads.  Usually, selecting the number of threads equal to the number of
cores (not HyperThreads) will produce the optimal performance.
