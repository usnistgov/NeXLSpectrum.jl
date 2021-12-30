# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis
```@meta
CurrentModule = NeXLSpectrum
```

# Spectrum Manipulation
See the [Spectrum Methods](@ref spectrum_methods) page for the most used methods and details.
```@docs
NeXLSpectrum.property!
NeXLSpectrum.duane_hunt
NeXLSpectrum.sigma
NeXLSpectrum.findsimilar
NeXLSpectrum.multiscore
NeXLSpectrum.multirank 
NeXLSpectrum.plot_compare
NeXLSpectrum.apply
NeXLSpectrum.loadmultispec
NeXLUncertainties.uv
NeXLSpectrum.χ²
NeXLSpectrum.recalibrate
NeXLSpectrum.shift
NeXLSpectrum.offset
NeXLSpectrum.dosenormalize
NeXLSpectrum.extent
NeXLSpectrum.characteristiccounts
NeXLSpectrum.scale
NeXLSpectrum.channelcount
NeXLSpectrum.sumcounts
NeXLSpectrum.shannon_entropy(::Spectrum)
```

## Spectrum Plotting
```@docs
Gadfly.plot( ::AbstractVector{Spectrum}; klms, edges, escapes, coincidences, autoklms, xmin, xmax, norm, yscale, ytransform, style, palette)
```

These types define the different ways that spectra can be scaled when plotted
using the `Gadfly.plot(...)` methods.
```@docs
NeXLSpectrum.SpectrumScaling
NeXLSpectrum.NoScaling
NeXLSpectrum.ScaleSum
NeXLSpectrum.ScaleDose
NeXLSpectrum.ScaleDoseWidth
NeXLSpectrum.ScaleROISum
NeXLSpectrum.ScalePeak
NeXLSpectrum.ScaleWidth
```

```@docs
NeXLSpectrum.NeXLSpectrumStyle
```

## Spectrum Tabulation
```@docs
NeXLUncertainties.asa(::Type{DataFrame}, ::Spectrum; properties)
NeXLUncertainties.asa(::Type{DataFrame}, ::AbstractArray{<:Spectrum})
```

# HyperSpectrum Manipulation
```@docs
NeXLSpectrum.HyperSpectrum
NeXLSpectrum.linescan
NeXLSpectrum.block
NeXLSpectrum.readrplraw
NeXLSpectrum.readptx
NeXLSpectrum.readhspy
NeXLSpectrum.ishspy
NeXLSpectrum.plane
NeXLSpectrum.roiimage
NeXLSpectrum.compress
NeXLSpectrum.maxpixel
NeXLSpectrum.colorize
NeXLSpectrum.labeledimages
NeXLSpectrum.labeledimage
NeXLSpectrum.region
NeXLSpectrum.indexofmaxpixel
NeXLSpectrum.roiimages
NeXLSpectrum.livetime!
```

# Fitting Filter
Core methods for constructing `FilterFitPacket`s and fitting spectra.
```@docs
NeXLSpectrum.reference
NeXLSpectrum.references
NeXLSpectrum.FilterFitPacket
NeXLSpectrum.suitability
NeXLSpectrum.suitablefor
NeXLSpectrum.fit_spectrum
NeXLSpectrum.FilteredReference
NeXLSpectrum.spectra
NeXLCore.elms
NeXLSpectrum.missingReferences
NeXLSpectrum.BasicFitResult
NeXLSpectrum.FilterFitResult
NeXLSpectrum.kratios
NeXLSpectrum.spectrum
NeXLSpectrum.residual
NeXLUncertainties.covariance
NeXLSpectrum.filteredresidual
NeXLUncertainties.extract
```

## Filter Fit Tabulation
```@docs
NeXLUncertainties.asa(::Type{DataFrame}, ::FilterFitPacket)
NeXLUncertainties.asa(::Type{DataFrame}, ::AbstractVector{<:FitResult}; charOnly, withUnc, format)
NeXLUncertainties.asa(::Type{DataFrame}, ::FilterFitResult; charOnly, material, columns, mc, fc)
NeXLUncertainties.asa(::Type{DataFrame}, ::FitResult; withUnc)
```

## Filter Fit Plotting
```@docs
Gadfly.plot(vq::VectorQuant, chs::UnitRange)
Gadfly.plot(ff::TopHatFilter, fr::FilteredReference)
Gadfly.plot(fr::FilteredReference; palette) 
Gadfly.plot(ffp::FilterFitPacket; kwargs...)
Gadfly.plot(ffr::FilterFitResult, roi::Union{Nothing, AbstractUnitRange{<:Integer}}; palette, style, xmax, comp, det, resp, yscale)
```

## Advanced Filter Fitting
```@docs
NeXLSpectrum.TopHatFilter
NeXLSpectrum.ConstantWidthFilter
NeXLSpectrum.GaussianFilter
NeXLSpectrum.VariableWidthFilter
NeXLSpectrum.tophatfilter
NeXLSpectrum.buildfilter
NeXLSpectrum.FilteredUnknownW
NeXLSpectrum.filterfit
NeXLSpectrum.isvisible
NeXLSpectrum.ReferenceLabel
NeXLSpectrum.SpectrumFeature
NeXLSpectrum.CharXRayLabel
NeXLSpectrum.EscapeLabel
NeXLSpectrum.UnknownLabel
NeXLSpectrum.charXRayLabels
```

# Matrix Correction
```@docs
NeXLMatrixCorrection.quantify
NeXLMatrixCorrection.estimatecoating
```

## Standardization
```@docs
NeXLSpectrum.extractStandards
NeXLCore.standardize
```

# EDS Detectors
```@docs
NeXLSpectrum.Detector
NeXLSpectrum.EDSDetector
NeXLSpectrum.resolution
NeXLSpectrum.Resolution
NeXLSpectrum.simpleEDSwICC
NeXLSpectrum.MnKaResolution
NeXLSpectrum.Beryllium
NeXLSpectrum.AP33Tabulation
NeXLSpectrum.matches
```

# EDS Detector Plotting
```@docs
Gadfly.plot(deteff::DetectorEfficiency)
Gadfly.plot(deteff::DetectorEfficiency, emax)
```

## Energy Axis Scales for EDS Detectors
Structures and functions that implement the energy scale functions for EDS detectors.
```@docs
NeXLSpectrum.EnergyScale
NeXLSpectrum.LinearEnergyScale
NeXLSpectrum.PolyEnergyScale
NeXLSpectrum.fwhm
NeXLSpectrum.gaussianwidth
NeXLSpectrum.resolution_to_fwhm
```

## Multi-Detector Functions
Functions for interpreting and manipulating spectr collected on multiple detectors simultaneously.
```@docs
NeXLSpectrum.multicompare
NeXLSpectrum.multimean
NeXLSpectrum.multisum
NeXLSpectrum.plot_multicompare
```

## Bremsstrahlung
```@docs
NeXLSpectrum.continuumrois
NeXLSpectrum.generated
NeXLSpectrum.emitted
NeXLSpectrum.fitcontinuum
NeXLSpectrum.detectorresponse
NeXLCore.weight
NeXLSpectrum.extents
NeXLSpectrum.profile
NeXLSpectrum.subtractcontinuum
NeXLSpectrum.heterogeneity
```

## Utility
```@docs
NeXLSpectrum.drawline
NeXLSpectrum.shannon_entropy(::AbstractArray{ColorTypes.Gray{FixedPointNumbers.N0f8}})
NeXLSpectrum.readSEManticsImage
NeXLCore.requiredbutmissing
NeXLCore.hasminrequired
NeXLSpectrum.annotate
```