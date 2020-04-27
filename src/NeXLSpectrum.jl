module NeXLSpectrum

using Reexport
using Requires

@reexport using NeXLCore

include("detector.jl")
export EnergyScale # Abstract: Detector energy calibration function
export LinearEnergyScale # A linear detector energy calibration function
export Resolution # Abstract: Detector resolution function
export MnKaResolution # Fiori function resolution function
export Detector # Abstract: X-ray detector
export EscapeArtifact # An escape peak
export SpectrumFeature # A CharXRay or EscapeArtifact (ComptonArtifact)
export EDSDetector # Abstract base for SimpleEDS and BasicEDS
export SimpleEDS # A simple EDS detector model
export BasicEDS # Adds min visibility by line family to SimpleEDS
export channel # energy to channel
export visible # Is a SpectrumFeature visible in the spectrum??
#export NeXLCore.energy # channel to energy
export linewidth # energy to linewidth
export channelcount # Detector channel count
export scale # Detector EnergyScale
export resolution # Detector Resolution
export simpleEDS # Create a basic EDS detector
export simpleEDSwICC # Create a basic EDS detector with incomplete charge collection
export extent # Determine the energy extent of x-ray lines on a detector ( Emin, Emax )
export extents # Determines contiguous channel extents from a set of characteristic lines or an element
export labeledextents # Like extents but labeled with a vector of the characteristic x-ray lines in each extent
export matching # Build a detector to match a spectrum
export lld # Low-level discriminator in eV

# Items defined in NeXL/spectrum.jl
include("spectrum.jl")
export Spectrum
#export NeXLCore.name # A human friendly name for the spectrum
export channel # channel for energy
export rangeofenergies # range of energies for channel `ch`
export channelwidth # width of channel ch
#export NeXLCore.energy # energy for channel
export dose # Spectrum probe dose
export counts # Spectrum channel data
export integrate # Integrate range of channels
export energyscale # Energy scale function
export subsample # Sub-sample a spectrum
export subdivide # Divide the counts among n spectra
export modelbackground # Model a background region
export modelBackground
export extractcharacteristic # Extract the characteristic intensity
export details # Outputs useful details about a spectrum
export peak # Estimates the peak intensity
export background # Estimates the background intensity
export estkratio # Estimate the k-ratio from two spectra for a ROI
export normalizedosewidth # Normalize intensity data to 1 nA⋅s⋅eV
export commonproperties
export maxspectrum

export maxproperty, minproperty # Min value of a property over a vector of spectra
export sameproperty # Returns the property value if all the spectra share the same value, errors otherwise

export minproperties # A list of the minimum required properties
export hasminrequired # Checks whether a spectrum has necessary properties
export requiredbutmissing # Lists missing properties

include("emsa.jl")
export readEMSA # Read an EMSA file
export writeEMSA # Write an EMSA file
# Also implements FileIO save(...) and load(...) for "ISO EMSA"
include("aspextiff.jl")
export readAspexTIFF
# Also implements FileIO save(...) and load(...) for "ASPEX TIFF"

include("hyperspectrum.jl")
export Signal  # The base class that makes hyperspectral data look like an Array of Real
export HyperSpectrum # The wrapper that makes a Signal look like an Array of Spectrum
export ashyperspectrum # Converts a Signal into a HyperSpectrum
export plane # Sum planes in a HyperSpectrum
export countmap # Convert a range of data channels into a Gray-scale image
export countmaps # Convert a vector of ranges-of-channels into Gray-scale images
export compressed # Compresses Integer type data down to the smallest size that will hold the max value.
export maxpixel # Bright's max-pixel derived spectrum
export indexofmaxpixel # Index producing the max pixel

include("rplraw.jl")
export RPLHeader
export readrplraw # Read a RPL/RAW file into a Signal
export readrpl # Read the RPL file into a RPLHeader
export writerplraw # Write a Signal into a RPL/RAW

include("fitlabels.jl")
export UnknownLabel#(spec)
export CharXRayLabel#(spec,roi,xrays)
export EscapeLabel#(spec,roi,xrays)
export charXRayLabels # Constructs CharXRayLabel(s)

include("filter.jl")
export TopHatFilter # Struct representing a fitting filter
export VariableWidthFilter # The default filter definition that varies the filter width with x-ray energy
export ConstantWidthFilter # An alternative filter definition that holds the filter width constant
export GaussianFilter # An alternative filter definition based on an offset Gaussian

export FilteredReference # A filtered datum representing a contiguous region of filtered reference data

# Different types of spectrum features that can be fit ROIs
export charFeature    # Characteristic peaks
export comptonFeature # Compton shifted characteristic peaks
export escapeFeature  # Escape artifact peaks

export buildfilter
export estimatebackground
export ascontiguous
export fitcontiguousg, fitcontiguousp, fitcontiguousw, fitcontiguouso
export filterfit
export filteredresidual

include("fitresult.jl")
export FilterFitResult
export findlabel
export fit
export heterogeneity
export peaktobackground
export characteristiccounts
export residual
export spectrum
export unknown
export kratios
export kratio

include("llsq.jl")
# The implementation for weighted filter-fit.
include("filterfit_wls.jl")
export FilteredUnknownW # A filtered datum representing an unknown spectrum (for weighted least squares fitting)
# The implementation for generalized filter-fit.
include("filterfit_gls.jl")
export FilteredUnknownG # A filtered datum representing an unknown spectrum (for generalized least squares fitting)

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
end

end
