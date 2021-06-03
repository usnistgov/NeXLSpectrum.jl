# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis

# Spectrum
```@docs
NeXLSpectrum.property!
NeXLSpectrum.duane_hunt
NeXLSpectrum.sigma
NeXLSpectrum.findsimilar
NeXLSpectrum.plot_compare  
NeXLSpectrum.multiscore
NeXLSpectrum.multirank 
NeXLSpectrum.apply
NeXLSpectrum.loadmultispec
NeXLUncertainties.uv
NeXLSpectrum.χ²
```

# Fitting Filter
```@docs
NeXLSpectrum.TopHatFilter
NeXLSpectrum.ConstantWidthFilter
NeXLSpectrum.GaussianFilter
NeXLSpectrum.VariableWidthFilter
NeXLSpectrum.tophatfilter
NeXLSpectrum.buildfilter
```

```@docs
NeXLSpectrum.FilteredReference
NeXLSpectrum.FilterFitPacket
NeXLSpectrum.spectra
NeXLCore.elms
NeXLSpectrum.reference
NeXLSpectrum.references
```

```@docs
NeXLSpectrum.FilteredUnknownW
NeXLSpectrum.FilteredUnknownG
NeXLSpectrum.filterfit
NeXLSpectrum.isvisible
NeXLSpectrum.fit_spectrum
NeXLSpectrum.missingReferences
```

```@docs
NeXLSpectrum.BasicFitResult
NeXLSpectrum.FilterFitResult
NeXLSpectrum.kratios
NeXLSpectrum.spectrum
NeXLSpectrum.residual
NeXLUncertainties.covariance
NeXLSpectrum.filteredresidual
```

```@docs
NeXLSpectrum.ReferenceLabel
NeXLSpectrum.SpectrumFeature
NeXLSpectrum.CharXRayLabel
NeXLSpectrum.EscapeLabel
NeXLSpectrum.UnknownLabel
```

## Matrix Correction
```@docs
NeXLMatrixCorrection.quantify
NeXLMatrixCorrection.estimatecoating
```

## Standardization
```@docs
NeXLSpectrum.suitability
NeXLSpectrum.suitablefor
NeXLSpectrum.extractStandards
NeXLCore.standardize
```

## Plot Scaling Modes
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
## Energy Axis Scales
```@docs
NeXLSpectrum.EnergyScale
NeXLSpectrum.LinearEnergyScale
NeXLSpectrum.PolyEnergyScale
```

```@docs
NeXLSpectrum.Detector
NeXLSpectrum.EDSDetector

NeXLSpectrum.resolution
NeXLSpectrum.Resolution
NeXLSpectrum.simpleEDSwICC
```

```@docs
NeXLSpectrum.HyperSpectrum
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

NeXLSpectrum.generated
NeXLSpectrum.continuumrois
NeXLSpectrum.emitted
NeXLSpectrum.fitcontinuum
NeXLSpectrum.indexofmaxpixel
NeXLSpectrum.roiimages
NeXLSpectrum.detectorresponse
NeXLCore.weight
NeXLSpectrum.extents
NeXLSpectrum.profile
NeXLSpectrum.region
NeXLSpectrum.livetime!

NeXLSpectrum.subtractcontinuum
NeXLSpectrum.heterogeneity
NeXLUncertainties.extract
NeXLSpectrum.extent
NeXLSpectrum.characteristiccounts
NeXLSpectrum.scale
NeXLSpectrum.charXRayLabels
NeXLSpectrum.channelcount
NeXLSpectrum.sumcounts

```

```@docs
NeXLSpectrum.Beryllium
NeXLSpectrum.AP33Tabulation
```

## Utility
```@docs
NeXLSpectrum.drawline
NeXLSpectrum.matches
```
