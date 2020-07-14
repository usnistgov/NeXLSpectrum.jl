# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis

`NeXLSpectrum` is a library of tools for manipulating EDS spectrum within the
NeXL framework. `NeXLSpectrum` depends on `NeXLUncertainties`, `NeXLCore` and
`NeXLMatrixCorrection` and loading `NeXLSpectrum` will also make these
libraries available.

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
