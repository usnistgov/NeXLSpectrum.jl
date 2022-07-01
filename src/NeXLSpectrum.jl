module NeXLSpectrum

using ThreadsX: maximum
using NeXLCore: isinside
using Base: UInt16, UInt64
using Reexport
using Requires

using Dates
using Images
using ImageAxes: AxisArrays, AxisArray, axes, axisnames, axisvalues
using Unitful: mm
using FileIO
using Mmap
using EzXML: readxml
using Polynomials: ImmutablePolynomial, fit, printpoly, roots, derivative, coeffs, Polynomials
using LinearAlgebra
using LoopVectorization
using Statistics
using DataAPI
using CSV
using DataFrames
using Interpolations: LinearInterpolation, AbstractInterpolation, bounds
using LsqFit: curve_fit
using ThreadsX

@reexport using NeXLCore
@reexport using NeXLMatrixCorrection

include("features.jl")
include("detector.jl")
export EnergyScale # Abstract: Detector energy calibration function
export LinearEnergyScale # A linear detector energy calibration function
export PolyEnergyScale # A polynomial energy calibration function
export Resolution # Abstract: Detector resolution function
export MnKaResolution # Fiori function resolution function
export Detector # Abstract: X-ray detector
export EscapeArtifact # An escape peak
export SpectrumFeature # A CharXRay or EscapeArtifact (ComptonArtifact)
export EDSDetector # Abstract base for BasicEDS and other EDS detector-like things
export BasicEDS # A basic EDS detector with min visibility by line family
export channel # energy to channel
export isvisible # Is a SpectrumFeature visible in the spectrum??
export linewidth # energy to linewidth
export channelcount # Detector channel count
export apply # Create a copy of a spectrum with the specified detector
export scale # Detector EnergyScale
export resolution # Detector Resolution
export resolution_to_fwhm # Given resolution at E what is resolution at Mn Kα
export profile # Computes the resolution function
export simpleEDS # Create a basic EDS detector
export simpleEDSwICC # Create a basic EDS detector with a naive incomplete charge collection model
export extent # Determine the energy extent of an x-ray line on a detector ( Emin, Emax )
export extents # Determine the energy extents of many x-ray lines on a detector
export matching # Build a detector to match a spectrum
export matches # Do the spectrum and detector match
export lld # Low-level discriminator in eV
export gaussianwidth, fwhm


# X-ray window models
include("window.jl")
# primary function transmission(wnd, energy, angle)
export AbstractWindow, ModeledWindow, TabulatedWindow, NoWindow
export WindowType, MoxtekAP33, MoxtekAP5, AmetekC1, AmetekC2, AmptekC1, AmptekC2, BerylliumWindow

include("detefficiency.jl")
export DetectorEfficiency
export efficiency
export SDDEfficiency, SiLiEfficiency # Helpers to build DetectorEfficiency
export detectorresponse # Build a matrix that describes the detectors response to X-rays
export detect # Detect Dict(CharXRay=>Intensity)
export simulate # Simulate a spectrum as measured on a detector

# Items defined in NeXL/spectrum.jl
include("spectrum.jl")
export Spectrum
#export NeXLCore.name # A human friendly name for the spectrum
export channel # channel for energy
export rangeofenergies # range of energies for channel `ch`
export channelwidth # width of channel ch
#export NeXLCore.energy # energy for channel
export dose # Spectrum probe dose
export property! # Set a property (easily broadcastable...)
export counts # Spectrum channel data
export integrate # Integrate range of channels
export kratio # Naive peak integration estimate of the k-ratio
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
export commonproperties # Properties that a collection of spectra share in common
export maxspectrum # Dave Bright's max spectra derived spectrum
export suitablefor # Which ROIs is a set of elements suitable for as reference for the specified element?
export missingReferences # A Vector with missing ROIs in a FilterFitPacket
export dosenormalize # Rescale the channel data to correspond to a different electron dose.
export offset # Offset the counts data by Δcounts

export maxproperty, minproperty # Min value of a property over a vector of spectra
export sameproperty # Returns the property value if all the spectra share the same value, errors otherwise
export textplot # A quick way to visualize a spectrum

include("compare.jl")
export findsimilar # Find the spectra that are most similar to each other
export χ² # Compare spectra
export sigma # Computes the channel-by-channel dose corrected difference from the mean. 

include("recalibrate.jl")
export recalibrate # Change the LinearEnergyScale against which a spectrum is calibrated
export shift # Shift the channels in a spectrum by a specified energy amount.

include("duanehunt.jl")
export duane_hunt # Estimate the Duane-Hunt limit

include("hyperspectrum.jl")
export HyperSpectrum # The wrapper that makes an Array{<:Real,N} look like an Array{Spectrum,N-1}
export plane # Sum planes in a HyperSpectrum
export roiimage # Convert a range of data channels into a Gray-scale image
export roiimages # Convert a vector of ranges-of-channels into Gray-scale images
export compress # Compresses Integer type data down to the smallest integer type that that will hold the max value.
export maxpixel, minpixel # Bright's max-pixel derived spectrum
export indexofmaxpixel # Index producing the max pixel
export avgpixel # Average intensity in a pixel
export sumcounts # Array containing the total number of counts in at each pixel
export depth # Number of spectral or result planes
export region # Extract a region from within the HyperSpectrum as a HyperSpectrum
export properties # The HyperSpectrum properties (mutable)
export axisname # The name of the i-th axis
export axisvalue # The calibrated coordinate value for the pixel coordinate
export axisrange # range of coordinate values for the specified axis
export livetime, livetime! # Get/Set livetime on a per-pixel basis
export block # Reduce the size of a HyperSpectrum by summing together blocks of adjacent pixels
export linescan # Extract the data in a 2D HyperSpectrum along a line (with width) as a 1D HyperSpectrum. 

include("rplraw.jl")
export RPLHeader
export readrplraw # Read a RPL/RAW file into a Array
export readrpl # Read the RPL file into a RPLHeader
export writerplraw # Write a Array into a RPL/RAW

include("emsa.jl")
include("aspextiff.jl")
include("brukerpdz.jl")
include("brukerspx.jl")
include("fileiosupport.jl")
include("semanticsptx.jl")
include("hspy.jl")

export loadspectrum # Load a spectrum from IO or filename
export savespectrum # Save a spectrum to IO or filename to a format
export sniffspectrum # Determine spectrum file type
export isspectrumfile # Is a file a spectrum file?
export readptx # Read a SEMantics PTX file
export readhspy # Read a HyperSpy-style HDF5 file
export ishspy # Is the file a HyperSpy-style HDF5 file

# SpectrumFileType structs for `loadspectrum()` and `savespectrum()` support
export BrukerSPX # Read only, SEM/EDS & XRF format
export BrukerPDZ # Read only,  XRF format
export ISOEMSA   # Read / write, ISO/EMSA spectrum file format
export ASPEXTIFF # Read only, Vendor SEM/EDS format

include("fitlabels.jl")
export ReferenceLabel
export UnknownLabel#(spec) <: ReferenceLabel
export CharXRayLabel#(spec|props,roi,xrays) <: QuantifiableLabel
export EscapeLabel#(spec,roi,xrays) <: ReferenceLabel
export charXRayLabels # Constructs CharXRayLabel(s) <: ReferenceLabel

include("filter.jl")
export TopHatFilter # Struct representing a fitting filter
export VariableWidthFilter # The default filter definition that varies the filter width with x-ray energy
export ConstantWidthFilter # An alternative filter definition that holds the filter width constant
export GaussianFilter # An alternative filter definition based on an offset Gaussian
export tophatfilter # Apply a top-hat filter to produce a FilteredReference
export FilteredReference # A filtered datum representing a contiguous region of filtered reference data

export buildfilter
export estimatebackground
export ascontiguous
export filterfit
export filteredresidual
export filterreference
export filterreferences

include("fitresult.jl")
export FitResult
export BasicFitResult
export FilterFitResult
export findlabel
export fit_spectrum # Preferred for individual spectra
export fit_spectra # Preferred for multiple spectra
export heterogeneity
export peaktobackground
export characteristiccounts
export residual
export spectrum
export unknown
export kratios
export kratio
export extractStandards

include("standardize.jl")
# export NeXLCore.standardize

include("llsq.jl")
# The implementation for weighted filter-fit.
include("filterfit_wls.jl")
export FilteredUnknownW # A filtered datum representing an unknown spectrum (for weighted least squares fitting)

include("continuum.jl")
export ContinuumModel # Build a model to compute the continuum intensity
export emitted # Compute the emitted continuum emission
export generated # Compute the generated continuum emission
export fitcontinuum # Fit a continuum model to a spectrum with :BeamEnergy, :TakeOffAngle, :Composition properties
export subtractcontinuum # Automatically remove the continuum contribution from a Spectrum
export fittedcontinuum
export continuumrois # The range of channels associated with the continuum

# These primarily exist to scale spectra for presentation (plotting)
include("specscaling.jl")
export SpectrumScaling # Abstract base
export NoScaling # <: SpectrumScaling
export ScaleDoseWidth # <: SpectrumScaling
export ScaleDose # <: SpectrumScaling
export ScaleSum # <: SpectrumScaling
export ScaleROISum # <: SpectrumScaling
export ScalePeak # <: SpectrumScaling
export ScaleWidth # <: SpectrumScaling
export scaledcounts # Return spectrum channel data scaled according the SpectrumScaling model

include("reference.jl")
export FilterFitPacket
export reference
export references
export spectra
export suitability # Tabulates material suitability for use as a fitting reference

include("qquant.jl")
export VectorQuant

include("direct.jl")
export DirectReference
export DirectReferences
export DirectFitResult
export direct

include("labeled.jl")
export labeledimage # Displays an image and caption.
export labeledimages # Displays a grid of images and captions.
export plotandimage # Display an image the right of a Gadfly plot

include("line.jl")
export drawline
export drawray

include("quantify.jl")
export quantify

include("smoothing.jl")
export SavitzkyGolayFilter

include("multidet.jl") # Multi-detector support
export loadmultispec # Loads multiple related spectra.
export multiscore # Scores spectra relative to one another. 
export multirank # A single number value that scores the spectra
export multisum, multimean
export multicompare

export plot_multicompare # Plots a channel-by-channel ratio comparison
export plot_compare # Plots a statistical channel-by-channel comparison of spectra

include("semantics_image.jl")
export readSEManticsImage
export annotate
export shannon_entropy

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
    @require Weave = "44d3d7a6-8a23-5bf8-98c5-b353f8df5ec9" include("weavesupport.jl")
end

end
