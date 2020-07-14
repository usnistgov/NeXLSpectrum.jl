# Model and fit the continuum
using NeXLMatrixCorrection

struct ContinuumModel
  mat::Material
  e0::Float64
  takeoff::Float64
  mc::Type{<:MatrixCorrection}
  br::Type{<:NeXLBremsstrahlung}
  """
      ContinuumModel(mat::Material, e0::Float64, det::DetectorEfficiency, takeoff::Float64)

  Create a continuum model for the specified material, beam energy, detector and take-off angle.  Computes
  the *detected* quantity of continuum generated in the sample.
  """
  ContinuumModel(
    mat::Material,
    e0::Float64,
    takeoff::Float64;
    matrixcorrection::Type{<:MatrixCorrection} = Riveros1993,
    bremsstrahlung::Type{<:NeXLBremsstrahlung} = Castellano2004b,
  ) = new(mat, e0, takeoff, matrixcorrection, bremsstrahlung)
end

"""
    emitted(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum emitted from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
function emitted(cm::ContinuumModel, ea::Float64) #
  g = generated(cm, ea)
  return g > 0.0 ? g * correctcontinuum(continuumcorrection(cm.mc, cm.mat, ea, cm.e0), cm.takeoff) : 0.0
end

"""
    generated(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum generated from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
generated(cm::ContinuumModel, ea::Float64) = #
  ea <= cm.e0 ? bremsstrahlung(cm.br, ea, cm.e0, cm.mat) : 0.0


"""
    fitcontinuum(
      spec::Spectrum,
      resp::AbstractArray,
      rois::Vector{UnitRange};
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
      mc::Type{<:MatricCorrection} = Riveros1993,
    )

    Fit a continuum model to the specified range of channels (`rois`).  The `resp` argument is a matrix which describes
the detector response on a channel-by-channel basis.  It can be calculated from an `EDSDetector` and an
`DetectorEfficiency` using `resp = NeXLSpectrum.detectorresponse(det, eff)`.  The `Spectrum` object must have
the :Composition, :BeamEnergy and :TakeOffAngle properties defined.
"""
function fitcontinuum(
  spec::Spectrum,
  resp::AbstractArray{<:Real,2},
  rois::AbstractArray{<:UnitRange,1};
  brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
  mc::Type{<:MatrixCorrection} = Riveros1993,
)
  @assert haskey(spec, :Composition) "The fitcontinuum(...) function requires the spec[:Composition] property."
  @assert haskey(spec, :BeamEnergy) "The fitcontinuum(...) function requires the spec[:BeamEnergy] property."
  @assert haskey(spec, :TakeOffAngle) "The fitcontinuum(...) function requires the spec[:TakeOffAngle] property."
  cmod = ContinuumModel(
    spec[:Composition],
    spec[:BeamEnergy],
    spec[:TakeOffAngle],
    matrixcorrection = mc,
    bremsstrahlung = brem,
  )
  minEmod = max(50.0, energy(lld(spec), spec)) # Ensures
  model = resp * map(e -> e > minEmod ? emitted(cmod, e) : 0.0, energyscale(spec))
  k = sum(sum(model[roi]) for roi in rois) / sum(sum(counts(spec, roi, Float64)) for roi in rois)
  props = Dict{Symbol,Any}(
    :TakeOffAngle => spec[:TakeOffAngle],
    :BeamEnergy => spec[:BeamEnergy],
    :Composition => spec[:Composition],
    :K => k,
    :Name => "Brem[Global][$(spec[:Name])]",
  )
  return Spectrum(spec.energy, model / k, props)
end


"""
    continuumrois(elms, det::EDSDetector, minE::Float64, maxE::Float64)

Compute the ROIs for the contiguous continuum regions for the specified elements `elms` on an
`EDSDetector` for the specified range of energies.
"""
function continuumrois(elms, det::EDSDetector, minE::Float64, maxE::Float64)::Vector{UnitRange{Int}}
  # Join contiguous rois
  function ascontiguous(rois)
    join(roi1, roi2) = min(roi1.start, roi2.start):max(roi1.stop, roi2.stop)
    srois = sort(rois)
    res = [srois[1]]
    for roi in srois[2:end]
      if length(intersect(res[end], roi)) > 0
        res[end] = join(roi, res[end])
      else
        push!(res, roi)
      end
    end
    return res
  end
  # Calculate the un-rois
  function invert(rois, over)
    st, res = over.start, []
    for roi in rois
      low, hgh = max(st, over.start), min(roi.start, over.stop)
      if hgh > low
        push!(res, low:hgh)
      end
      st = roi.stop
    end
    if st < over.stop
      push!(res, st:over.stop)
    end
    return res
  end
  extend(roi, n) = roi.start-n:roi.stop+n
  # Average channel width between minE and maxE.
  function avgwidth(det, minE, maxE)
    minC, maxC = channel(minE,det), channel(maxE,det)
    return (energy(maxC,det) - energy(minC,det))/(maxC-minC)
  end
  extra = round(Int, (2.0 * resolution(enx"Mn K-L3", det)) / avgwidth(det,minE,maxE))
  rois = mapreduce(elm -> extend.(extents(elm, det, 1.0e-3), extra), append!, elms)
  return invert(ascontiguous(rois), channel(minE, det):channel(maxE, det))
end



"""
    fitcontinuum(
      spec::Spectrum,
      det::EDSDetector,
      resp::AbstractArray{<:Real,2}; #
      minE::Float64 = 1.5e3,
      maxE::Float64 = 0.95 * spec[:BeamEnergy],
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
      mc::Type{<:MatrixCorrection} = Riveros1993,
    )

Fit the continuum from ROIs determined from the data within the spectrum (:Composition, :BeamEnergy & :TakeOffAngle).
The ROIs are computed using `continuumrois(...)` and each roi is fit seperately.


    fittedcontinuum(
      spec::Spectrum,
      det::EDSDetector,
      resp::AbstractArray{<:Real,2}; #
      mode = :Global [ | :Local ] # Fit to all ROIs simultaneously (:Global) or to each roi independently (:Local)
      minE::Float64 = 1.5e3,
      maxE::Float64 = 0.95 * spec[:BeamEnergy],
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
      mc::Type{<:MatrixCorrection} = Riveros1993,
    )::Spectrum

Fit the continuum under the characteristic peaks by fitting the closest continuum ROIs.  The low energy peaks are
fit using the continuum immediately higher in energy and the high energy peaks are fit using the continuum on both
sides.
"""
function fitcontinuum(
  spec::Spectrum,
  det::EDSDetector,
  resp::AbstractArray{<:Real,2}; #
  minE::Float64 = 1.5e3,
  maxE::Float64 = 0.95 * spec[:BeamEnergy],
  brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
  mc::Type{<:MatrixCorrection} = Riveros1993,
)::Vector{Float64}
  @assert haskey(spec, :Composition) "The fitcontinuum(...) function requires the spec[:Composition] property."
  @assert haskey(spec, :BeamEnergy) "The fitcontinuum(...) function requires the spec[:BeamEnergy] property."
  @assert haskey(spec, :TakeOffAngle) "The fitcontinuum(...) function requires the spec[:TakeOffAngle] property."
  cmod = ContinuumModel(
    spec[:Composition],
    spec[:BeamEnergy],
    spec[:TakeOffAngle],
    matrixcorrection = mc,
    bremsstrahlung = brem,
  )
  # meas is the raw continuum shape.  It needs to be scaled to the unknown.
  minEmod = max(50.0, energy(lld(spec), spec))
  model = resp * map(e -> e > minEmod ? emitted(cmod, e) : 0.0, energyscale(spec))
  prevroi, brem = missing, counts(spec, Float64)
  for roi in continuumrois(elms(spec), det, minE, maxE)
    currrois = ismissing(prevroi) ? (roi,) : (prevroi, roi)
    k = sum(sum(counts(spec, roi, Float64)) for roi in currrois) / sum(sum(model[roi]) for roi in currrois)
    peak = (ismissing(prevroi) ? lld(det) : prevroi.stop):roi.start
    brem[peak] = k * model[peak]
    prevroi = roi
  end
  return brem
end

function fittedcontinuum(
  spec::Spectrum,
  det::EDSDetector,
  resp::AbstractArray{<:Real,2}; #
  mode::Symbol = :Global,
  minE::Float64 = 1.5e3,
  maxE::Float64 = 0.95 * spec[:BeamEnergy],
  brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
  mc::Type{<:MatrixCorrection} = Riveros1993,
)::Spectrum
  function localfittedcontinuum(spec, det, resp, minE, maxE, brem, mc)::Spectrum
    res = Spectrum(
      spec.energy,
      fitcontinuum(spec, det, resp, minE = minE, maxE = maxE, brem = brem, mc = mc),
      copy(spec.properties),
    )
    res[:Name] = "Brem[Local][$(spec[:Name])]"
    return res
  end
  globalfittedcontinuum(spec, det, resp, minE, maxE, brem, mc)::Spectrum =
    fitcontinuum(spec, resp, continuumrois(elms(spec), det, minE, maxE), brem=brem, mc=mc)
  @assert (mode==:Global) || (mode==:Local) "mode must equal :Global | :Local in fitted continuum"
  return mode == :Global ? globalfittedcontinuum(spec, det, resp, minE, maxE, brem, mc) : #
    localfittedcontinuum(spec, det, resp, minE, maxE, brem, mc)
end

"""
    subtractcontinuum(
      spec::Spectrum,
      det::EDSDetector,
      resp::AbstractArray{<:Real,2}; #
      minE::Float64 = 1.5e3,
      maxE::Float64 = 0.95 * spec[:BeamEnergy],
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
      mc::Type{<:MatrixCorrection} = Riveros1993,
    )::Spectrum

Computes the characteristic-only spectrum by subtracting the .
"""
function subtractcontinuum(
  spec::Spectrum,
  det::EDSDetector,
  resp::AbstractArray{<:Real,2}; #
  minE::Float64 = 1.5e3,
  maxE::Float64 = 0.95 * spec[:BeamEnergy],
  brem::Type{<:NeXLBremsstrahlung} = Castellano2004a,
  mc::Type{<:MatrixCorrection} = Riveros1993,
)::Spectrum
  res = Spectrum(
    spec.energy,
    counts(spec, Float64) - fitcontinuum(spec, det, resp, minE = minE, maxE = maxE, brem = brem, mc = mc),
    copy(spec.properties),
  )
  res[:Name] = "CharOnly[$(res[:Name])]"
  return res
end
