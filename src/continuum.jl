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
  ContinuumModel(mat::Material, e0::Float64, det::DetectorEfficiency, takeoff::Float64) =
      new(mat, e0, takeoff, det, CitZAF, Castellano2004a)
end

"""
    emitted(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum emitted from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
emitted(cm::ContinuumModel, ea::Float64) = #
  generated(cm, ea) * correctcontinuum(continuumcorrection(cm.mc, cm.mat, ea, cm.e0), ea, cm.takeoff)


"""
    generated(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum generated from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
generated(cm::ContinuumModel, ea::Float64) = #
  bremsstrahlung(cm.br, ea, cm.e0, cm.mat) * efficiency(cm.detector, ea, π/2)
