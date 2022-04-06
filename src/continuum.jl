# Model and fit the continuum

struct ContinuumModel
    material::Material
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
        material::Material,
        e0::Float64,
        takeoff::Float64;
        matrixcorrection::Type{<:MatrixCorrection} = Riveros1993,
        bremsstrahlung::Type{<:NeXLBremsstrahlung} = Castellano2004b,
    ) = new(material, e0, takeoff, matrixcorrection, bremsstrahlung)
end

"""
    emitted(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum emitted from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
function emitted(cm::ContinuumModel, ea::Float64) #
    g = generated(cm, ea)
    return g > 0.0 ?
           g * correctcontinuum(
        continuumcorrection(cm.mc, cm.material, ea, cm.e0),
        cm.takeoff,
    ) : 0.0
end

"""
    generated(cm::ContinuumModel, ea::Float64)

Compute the intensity of the measured continuum generated from the material and conditions specified in the continuum
model object at the specified measured energy `ea`.
"""
generated(cm::ContinuumModel, ea::Float64) = #
    ea <= cm.e0 ? bremsstrahlung(cm.br, ea, cm.e0, cm.material) : 0.0


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
    k =
        sum(sum(model[roi]) for roi in rois) /
        sum(sum(counts(spec, roi, Float64)) for roi in rois)
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
    continuumrois(elems, det::EDSDetector, minE::Float64, maxE::Float64)

Compute the ROIs for the contiguous continuum regions for the specified elements `elems` on an
`EDSDetector` for the specified range of energies.
"""
function continuumrois(
    elems,
    det::EDSDetector,
    minE::Float64,
    maxE::Float64,
)::Vector{UnitRange{Int}}
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
    avgwidth = let 
        minC, maxC = channel(minE, det), channel(maxE, det)
        (energy(maxC, det) - energy(minC, det)) / (maxC - minC)
    end
    extra = round(Int, (2.0 * resolution(enx"Mn K-L3", det)) / avgwidth)
    rois = mapreduce(elm -> extend.(extents(elm, det, 1.0e-3), extra), append!, elems)
    # Join rois into contiguous rois
    ascontiguous = let
        srois = sort(rois)
        res = [srois[1]]
        for roi in srois[2:end]
            if length(intersect(res[end], roi)) > 0
                res[end] = min(roi.start, res[end].start):max(roi.stop, res[end].stop)
            else
                push!(res, roi)
            end
        end
        res
    end
    return invert(ascontiguous, channel(minE, det):channel(maxE, det))
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
      width::Int = 20, # Width of ROI at each end of each patch of continuum that is matched
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
    maxE::Float64 = 0.90 * spec[:BeamEnergy],
    width::Int = 20,
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
    chlow, k0 = 1, 0.0
    chdata, result = counts(spec, Float64), zeros(Float64, length(spec))
    for roi in continuumrois(elms(spec, true), det, minE, maxE)
        # Begining of continuum ROI
        roi1 = roi.start:(roi.start+min(width, length(roi))-1)
        k1 = sum(@view chdata[roi1]) / sum(@view model[roi1])
        # End of continuum ROI
        roi2 = (roi.stop - min(width, length(roi)) + 1):roi.stop
        k2 = sum(@view chdata[roi2]) / sum(@view model[roi2])
        # Model between roi1 and roi2 
        k1_2(ch) = (k2-k1)*(ch-roi.start)/length(roi) + k1
        result[roi] = k1_2.(roi) .* model[roi]
        k0 = chlow == 1 ? k1 : k0
        # Model between prev roi and roi
        roi0_1 = chlow:(roi.start-1)
        k0_1(ch) = (k1-k0)*(ch-roi0_1.start)/length(roi0_1) + k0
        result[roi0_1] = k0_1.(roi0_1) .* model[roi0_1]
        chlow, k0 = roi.stop + 1, k2
    end
    # Fill the final patch
    result[chlow:length(result)] = k0*model[chlow:length(result)]
    return result
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
    globalfittedcontinuum(spec, det, resp, minE, maxE, brem, mc)::Spectrum = fitcontinuum(
        spec,
        resp,
        continuumrois(elms(spec, true), det, minE, maxE),
        brem = brem,
        mc = mc,
    )
    @assert (mode == :Global) || (mode == :Local) "mode must equal :Global | :Local in fitted continuum"
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

Computes the characteristic-only spectrum by subtracting the continuum.
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
        counts(spec, Float64) -
        fitcontinuum(spec, det, resp, minE = minE, maxE = maxE, brem = brem, mc = mc),
        copy(spec.properties),
    )
    res[:Name] = "CharOnly[$(res[:Name])]"
    return res
end
