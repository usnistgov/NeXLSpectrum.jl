# ![](NeXL_sm.png)Spectrum
## Microanalytical X-ray Spectrum Analysis

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
```

```@docs
NeXLSpectrum.FilteredUnknownW
NeXLSpectrum.FilteredUnknownG
NeXLSpectrum.filterfit

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

NeXLSpectrum.plane
NeXLSpectrum.roiimage
NeXLSpectrum.asimage
NeXLSpectrum.compressed
NeXLSpectrum.maxpixel

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


NeXLSpectrum.subtractcontinuum
NeXLSpectrum.heterogeneity
NeXLSpectrum.visible
NeXLUncertainties.extract
NeXLSpectrum.extent
NeXLSpectrum.characteristiccounts
NeXLSpectrum.scale
NeXLSpectrum.HyperspectrumQuant
NeXLSpectrum.charXRayLabels
NeXLSpectrum.channelcount
```

```@docs
NeXLSpectrum.Beryllium
NeXLSpectrum.AP33Tabulation
```
