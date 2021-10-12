# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis

# Spectrum
```@docs
NeXLSpectrum.Spectrum
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
NeXLSpectrum.recalibrate
NeXLSpectrum.shift
NeXLSpectrum.offset
NeXLSpectrum.dosenormalize
NeXLSpectrum.extent
NeXLSpectrum.characteristiccounts
NeXLSpectrum.scale
NeXLSpectrum.charXRayLabels
NeXLSpectrum.channelcount
NeXLSpectrum.sumcounts
```

# HyperSpectrum
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
```

# Matrix Correction
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
```
## Energy Axis Scales for EDS Detectors
```@docs
NeXLSpectrum.EnergyScale
NeXLSpectrum.LinearEnergyScale
NeXLSpectrum.PolyEnergyScale
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
NeXLSpectrum.matches
```