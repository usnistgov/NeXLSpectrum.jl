module NeXLSpectrum

using NeXLCore

include("detector.jl")
export EnergyScale # Abstract: Detector energy calibration function
export LinearEnergyScale # A linear detector energy calibration function
export Resolution # Abstract: Detector resolution function
export MnKaResolution # Fiori function resolution function
export Detector # Abstract: X-ray detector
export SimpleEDS # A simple EDS detector model
export channel # energy to channel
export energy # channel to energy
export linewidth # energy to linewidth
export channelcount # Detector channel count
export scale # Detector EnergyScale
export resolution # Detector Resolution
export basicEDS # Create a basic EDS detector
export extent # Determine the channel extent of lines and line families

# Items defined in NeXL/spectrum.jl
include("spectrum.jl")
export Spectrum
export name # A human friendly name for the spectrum
export channel # channel for energy
export energy # energy for channel
export readEMSA # Read an EMSA file
export dose # Spectrum probe dose
export counts # Spectrum channel data
export integrate # Integrate range of channels
export energyscale # Energy scale function
export subsample # Sub-sample a spectrum
export modelbackground # Model a background region
export modelBackground
export extractcharacteristic # Extract the characteristic intensityS

end
