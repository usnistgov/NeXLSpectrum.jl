using CubicSplines
# Model and fit the continuum

struct ContinuumModel
    material::Material
    e0::Float64
    takeoff::Float64
    mc::Type{<:MatrixCorrection}
    br::Type{<:NeXLBremsstrahlung}
    """
        ContinuumModel(
            material::Material,
            e0::Float64,
            takeoff::Float64;
            matrixcorrection::Type{<:MatrixCorrection}=CitZAF,
            bremsstrahlung::Type{<:NeXLBremsstrahlung}=Castellano2004b
        )

    Create a continuum model for the specified material, beam energy, detector and take-off angle.  Computes
    the *detected* quantity of continuum generated in the sample.
    """
    ContinuumModel(
        material::Material,
        e0::Float64,
        takeoff::Float64;
        matrixcorrection::Type{<:MatrixCorrection}=CitZAF,
        bremsstrahlung::Type{<:NeXLBremsstrahlung}=Castellano2004b
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
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004b,
      mc::Type{<:MatricCorrection} = CitZAF,
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
    brem::Type{<:NeXLBremsstrahlung}=Castellano2004b,
    mc::Type{<:MatrixCorrection}=CitZAF
)
    @assert haskey(spec, :Composition) "The fitcontinuum(...) function requires the spec[:Composition] property."
    @assert haskey(spec, :BeamEnergy) "The fitcontinuum(...) function requires the spec[:BeamEnergy] property."
    @assert haskey(spec, :TakeOffAngle) "The fitcontinuum(...) function requires the spec[:TakeOffAngle] property."
    cmod = ContinuumModel(
        spec[:Composition],
        spec[:BeamEnergy],
        spec[:TakeOffAngle],
        matrixcorrection=mc,
        bremsstrahlung=brem,
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
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004b,
      mc::Type{<:MatrixCorrection} = CitZAF,
    )

Fit the continuum from ROIs determined from the data within the spectrum (:Composition, :BeamEnergy & :TakeOffAngle).
The ROIs are computed using `continuumrois(...)` and each roi is fit seperately.
"""
function fitcontinuum(
    spec::Spectrum,
    det::EDSDetector,
    resp::AbstractArray{<:Real,2}; #
    minE::Float64=1.5e3,
    maxE::Float64=0.90 * spec[:BeamEnergy],
    width::Int=20,
    brem::Type{<:NeXLBremsstrahlung}=Castellano2004b,
    mc::Type{<:MatrixCorrection}=CitZAF
)
    @assert haskey(spec, :Composition) "The fitcontinuum(...) function requires the spec[:Composition] property."
    @assert haskey(spec, :BeamEnergy) "The fitcontinuum(...) function requires the spec[:BeamEnergy] property."
    @assert haskey(spec, :TakeOffAngle) "The fitcontinuum(...) function requires the spec[:TakeOffAngle] property."
    cmod = ContinuumModel(
        spec[:Composition],
        spec[:BeamEnergy],
        spec[:TakeOffAngle],
        matrixcorrection=mc,
        bremsstrahlung=brem,
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
        roi2 = (roi.stop-min(width, length(roi))+1):roi.stop
        k2 = sum(@view chdata[roi2]) / sum(@view model[roi2])
        # Model between roi1 and roi2 
        k1_2(ch) = (k2 - k1) * (ch - roi.start) / length(roi) + k1
        result[roi] = k1_2.(roi) .* model[roi]
        k0 = chlow == 1 ? k1 : k0
        # Model between prev roi and roi
        roi0_1 = chlow:(roi.start-1)
        k0_1(ch) = (k1 - k0) * (ch - roi0_1.start) / length(roi0_1) + k0
        result[roi0_1] = k0_1.(roi0_1) .* model[roi0_1]
        chlow, k0 = roi.stop + 1, k2
    end
    # Fill the final patch
    result[chlow:length(result)] = k0 * model[chlow:length(result)]
    props = Dict{Symbol,Any}(
        :TakeOffAngle => spec[:TakeOffAngle],
        :BeamEnergy => spec[:BeamEnergy],
        :Composition => spec[:Composition],
        :Name => "Brem[Local][$(spec[:Name])]",
    )
    return Spectrum(spec.energy, result, props)
end

"""
fittedcontinuum(
    spec::Spectrum,
    det::EDSDetector,
    resp::AbstractArray{<:Real,2}; #
    mode = :Global [ | :Local ] # Fit to all ROIs simultaneously (:Global) or to each roi independently (:Local)
    minE::Float64 = 1.5e3,
    maxE::Float64 = 0.95 * spec[:BeamEnergy],
    width::Int = 20, # Width of ROI at each end of each patch of continuum that is matched
    brem::Type{<:NeXLBremsstrahlung} = Castellano2004b,
    mc::Type{<:MatrixCorrection} = CitZAF,
    thresh::AbstractFloat=5.0, # Statistical threshold to be considered continuum
    width::Integer=10, # Additional width added to peak regions (in channels)
    filt=buildfilter(eltype(spec), VariableWidthFilter, det) # Filter applied to spectrum
  )::Spectrum

Fit the continuum under the characteristic peaks by fitting the closest continuum ROIs.  The low energy peaks are
fit using the continuum immediately higher in energy and the high energy peaks are fit using the continuum on both
sides.

  * mode = :Global [ | :Local ] Global fits the model to the data using a single scale factor. :Local tweaks the
  global model at ROIs above and below the characteristic peaks.

  Typically, :Global produces the overall best fit but :Local fits better around the characteristic peaks and is 
  better for modeling the continuum under the peaks.
"""
function fittedcontinuum(
    spec::Spectrum,
    det::EDSDetector,
    resp::AbstractArray{<:Real,2}; #
    mode::Symbol=:Global,
    minE::Float64 = 1.5e3,
    maxE::Float64=0.95 * spec[:BeamEnergy],
    brem::Type{<:NeXLBremsstrahlung}=Castellano2004b,
    mc::Type{<:MatrixCorrection}=CitZAF,
    thresh::AbstractFloat=5.0,
    width::Integer=10,
    filt=buildfilter(eltype(spec), VariableWidthFilter, det)
)
    r = 1:min(channel(maxE, det), length(spec))
    function thresholdrois(fs1, r; w=width, th=thresh)
        res = abs.(@view fs1.filtered[r]) .< th*sqrt.(@view fs1.covariance[r])
        max = r.start
        for i in r
            if !res[i]
                max = i + w
            end
            res[i] = (i > max)
        end
        min = r.stop
        for i in reverse(r)
            if !res[i]
                min = i - w
            end
            res[i] = (i < min)
        end
        rois = UnitRange{Int64}[]
        i=r.start
        while i < r.stop
            if res[i]
                j=i+1
                while (j<r.stop) && res[j]
                    j+=1
                end
                push!(rois, i:j)
                i=j
            end
            i+=1
        end
        return rois
    end
    @assert (mode == :Global) || (mode == :Local) "mode must equal :Global | :Local in fitted continuum"
    fs1=tophatfilter(spec,filt)
    crois = thresholdrois(fs1, r)
    gl = fitcontinuum(spec, resp, crois; brem=brem, mc=mc)
    return mode == :Global ? gl : tweakcontinuum(spec, gl, crois)
end

"""
    tweakcontinuum(
        meas::Spectrum{T}, 
        cont::Spectrum{U}, 
        crois::Vector{UnitRange{Int}}; 
        nch=10, 
        maxSc=1.5
    ) where { T<: Real, U <: Real }

Takes a measured spectrum and the associated modeled continuum spectrum and tweaks the continuum spectrum to produce
a better fit to the measured.  It focuses on the ends of the continuum ROIs (in `crois`) to match the continuum at
these channels.

  * `maxSc` limits how large a tweak is acceptable.
  * `nch` determines how many channels to match at the end of each ROI in `crois`
"""
function tweakcontinuum(meas::Spectrum{T}, cont::Spectrum{U}, crois::Vector{UnitRange{Int}}; nch=5, maxSc=2.0) where {T<:Real,U<:Real}
    function LinearSpline(x, y)
        @assert length(x)==length(y)
        function f(xx)
            i = findfirst(v -> xx <= v, x)
            return if !isnothing(i)
                i > firstindex(x) ? ((y[i] - y[i-1]) / (x[i] - x[i-1])) * (xx - x[i-1]) + y[i-1] : y[i]
            else
                y[lastindex(y)]
            end
        end
        return f
    end
    x, y, chmax = Float64[1.0], Float64[1.0], channel(0.9*meas[:BeamEnergy], meas)
    for croi in crois
        if (length(croi) > (3 * nch) / 2) && (croi.start < chmax)
            # Add a pivot beginning and end of croi
            stroi = croi.start:croi.start+nch-1
            push!(x, stroi.start)
            push!(y, Base.clamp(integrate(meas, stroi) / max(integrate(cont, stroi),1.0), 1.0 / maxSc, maxSc))
            if croi.stop < chmax
                endroi = croi.stop-nch+1:min(length(cont), croi.stop)
                push!(x, endroi.stop)
                push!(y, clamp(integrate(meas, endroi) / max(integrate(cont, endroi), 1.0), 1.0 / maxSc, maxSc))
            else
                push!(x, chmax)
                push!(y, 1.0)
            end
        end
    end
    if x[end] < length(cont)
        push!(x, length(cont))
        push!(y, 1.0)
    end
    # In case there are too few points in the spline for a cubic spline
    spline = length(x) >= 5 ? CubicSpline(x, y) : LinearSpline(x, y)
    scale = spline.(eachindex(cont.counts))
    props = Dict(
        cont.properties...,
        :CScale => scale,
        :Parent => cont,
        :Name => "Tweaked[$(cont[:Name])]"
    )
    return Spectrum(cont.energy, scale .* cont.counts, props)
end

"""
    subtractcontinuum(
      spec::Spectrum,
      det::EDSDetector,
      resp::AbstractArray{<:Real,2}; #
      minE::Float64 = 1.5e3,
      maxE::Float64 = 0.95 * spec[:BeamEnergy],
      brem::Type{<:NeXLBremsstrahlung} = Castellano2004b,
      mc::Type{<:MatrixCorrection} = CitZAF,
    )::Spectrum

Computes the characteristic-only spectrum by subtracting the continuum.
"""
function subtractcontinuum(
    spec::Spectrum,
    det::EDSDetector,
    resp::AbstractArray{<:Real,2}; #
    minE::Float64=1.5e3,
    maxE::Float64=0.95 * spec[:BeamEnergy],
    brem::Type{<:NeXLBremsstrahlung}=Castellano2004b,
    mc::Type{<:MatrixCorrection}=CitZAF
)::Spectrum
    res = Spectrum(
        spec.energy,
        counts(spec, Float64) -
        fitcontinuum(spec, det, resp, minE=minE, maxE=maxE, brem=brem, mc=mc),
        copy(spec.properties),
    )
    res[:Name] = "CharOnly[$(res[:Name])]"
    return res
end
