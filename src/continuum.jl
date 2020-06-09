# Model and fit the continuum
using NeXLMatrixCorrection

struct ContinuumModel
  mat::Material
  e0::Float64
  takeoff::Float64
  detector::DetectorEfficiency
  mc::Type{<:MatrixCorrection}
  br::Type{<:NeXLBremsstrahlung}
"""
    ContinuumModel(mat::Material, e0::Float64, det::DetectorEfficiency, takeoff::Float64)

Create a continuum model for the specified material, beam energy, detector and take-off angle.  Computes
the *detected* quantity of continuum generated in the sample.
"""
  ContinuumModel(mat::Material, e0::Float64, det::DetectorEfficiency, takeoff::Float64; matrixcorrection::Type{<:MatrixCorrection}=CitZAF, bremsstrahlung::Type{<:NeXLBremsstrahlung}=Castellano2004b) =
      new(mat, e0, takeoff, det, matrixcorrection, bremsstrahlung)
end

"""
    emitted(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum emitted from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
function emitted(cm::ContinuumModel, ea::Float64) #
  g = generated(cm, ea)
  return g > 0.0 ? g * correctcontinuum(continuumcorrection(cm.mc, cm.mat, ea, cm.e0), ea, cm.takeoff) : 0.0
end

"""
    generated(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum generated from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
generated(cm::ContinuumModel, ea::Float64) = #
  ea <= cm.e0 ? bremsstrahlung(cm.br, ea, cm.e0, cm.mat) * efficiency(cm.detector, ea, π/2) : 0.0


"""
    detectorresponse(det::EDSDetector, eff::DetectorEfficiency, incidence::Float64=π/2)::AbstractMatrix

Build a matrix which models the detector response including aspects like the detector efficiency, the resolution, the
escape peaks.  All the warts that can be modeled within a linear model but not things like coincidence peaks that
are non-linear.  This function can (!should!) be specialized for more sophisticated detector models that include more
warts.

Example:

    genint = computegeneratedintensity(....) # Either characteristic or Bremsstrahlung...
    det = simpleEDS(4096, 5.0, 0.0, 132.0, 10)
    eff = SDDEfficiency(AP33Tabulation(); thickness=0.0370, deadlayer=30.0e-7, entrance=Film(pure(n"Al"), 10.0e-7))
    resp = detectorresponse(det, eff)
    # finally compute the measured signal
    measured = genint*resp
"""
function detectorresponse(det::EDSDetector, eff::DetectorEfficiency, incidence::Float64=π/2)::AbstractMatrix
  res = zeros(Float64, (channelcount(det), channelcount(det)))
  # An x-ray with energy in ch will be dispersed among a range of channels about ch
  for ch in lld(det):channelcount(det) # X-ray energy by channel
      el, eh = energy(ch, det), energy(ch+1, det)  # X-ray energy
      effic, fwhm = 0.5*(efficiency(eff, eh, incidence) + efficiency(eff, el, incidence)), resolution(0.5*(eh+el), det)
      roc = max(lld(det), channel( el - 3.0 * fwhm, det)):min(channelcount(det), channel(el + 3.0 * fwhm, det))
      res[ch, roc] = map(ch2->effic * profile(ch2, 0.5*(eh+el), det), roc)
  end
  return res
end
