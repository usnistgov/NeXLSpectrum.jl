using Distributions: Poisson

# Model detector efficiency (Window +  Detector Crystal)

struct DetectorEfficiency
    name::String
    window::AbstractWindow
    surface::Vector{Film}
    active::Film
end

Base.show(io::IO, de::DetectorEfficiency) = print(io, "$(de.name)[$(de.window)]")

SDDEfficiency(
    window::AbstractWindow;
    thickness=0.0370,
    deadlayer=30.0e-7,
    entrance=Film(pure(n"Al"), 10.0e-7)
) = DetectorEfficiency(
    "SDD",
    window,
    [entrance, Film(pure(n"Si"), deadlayer)],
    Film(pure(n"Si"), thickness),
)

SiLiEfficiency(
    window::AbstractWindow;
    thickness=0.250,
    deadlayer=30.0e-7,
    entrance=Film(pure(n"Al"), 10.0e-7)
) = DetectorEfficiency(
    "Si(Li)",
    window,
    [entrance, Film(pure(n"Si"), deadlayer)],
    Film(pure(n"Si"), thickness),
)

function efficiency(aa::DetectorEfficiency, energy::Float64, angle::Float64=π / 2)
    return energy > 20.0 ? #
           transmission(aa.window, energy, angle) *
           (1.0 - transmission(aa.active, energy, angle)) *
           mapreduce(lyr -> transmission(lyr, energy, angle), *, aa.surface) : #
           0.0
end

"""
    detectorresponse(det::EDSDetector, eff::DetectorEfficiency, incidence::Float64=π/2)::AbstractMatrix

Build a matrix which models the detector response including aspects like the detector efficiency, the resolution, the
escape peaks.  All the warts that can be modeled within a linear model but not things like coincidence peaks that
are non-linear.  This function can (!should!) be specialized for more sophisticated detector models that include more
warts.

It also dicretizes the input energies on the same scale as the `EDSDetector` (thus it is square.)  This is reasonable
when the detector channel width is much less than the resolution.

Example:

    genint = computegeneratedintensity(....) # Either characteristic or Bremsstrahlung...
    det = simpleEDS(4096, 5.0, 0.0, 132.0, 10)
    eff = SDDEfficiency(AP33Tabulation(); thickness=0.0370, deadlayer=30.0e-7, entrance=Film(pure(n"Al"), 10.0e-7))
    resp = detectorresponse(det, eff)
    # finally compute the measured signal
    measured = resp*genint
"""
function detectorresponse(
    det::EDSDetector,
    eff::DetectorEfficiency,
    incidence::Float64=π / 2,
)
    res = zeros(Float64, (channelcount(det), channelcount(det)))
    # An x-ray with energy in ch will be dispersed among a range of channels about ch
    for ch in max(channel(10.0, det), lld(det)):channelcount(det) # X-ray energy by channel
        el, eh = energy(ch, det), energy(ch + 1, det)  # X-ray energy
        @assert el > 0.0
        fwhm = resolution(0.5 * (eh + el), det)
        @assert !isnan(fwhm)
        # Full width of detectable X-rays
        full_roc = channel(el - 3.0 * fwhm, det):channel(el + 3.0 * fwhm, det)
        # Range of available channels
        in_roc = max(lld(det), full_roc.start):min(channelcount(det), full_roc.stop)
        # Nominal efficiency for an energy(ch) X-ray
        effic = 0.5 * (efficiency(eff, eh, incidence) + efficiency(eff, el, incidence))
        @assert !isnan(effic)
        prof = map(ch2 -> profile(ch2, 0.5 * (eh + el), det), full_roc)
        pre = (in_roc.start - full_roc.start) + 1
        post = length(prof) - (full_roc.stop - in_roc.stop)
        res[in_roc, ch] = effic * @view prof[pre:post]
    end
    return res
end

"""
    detect(emitted::Dict{<:XRay,<:Real}, det::EDSDetector, response::Matrix{Float64}; noise=false)

Returns a Spectrum as though the intensities in `emitted` were detected on the specified detector.
`noise=true` will add Poisson noise to the resulting measured `Spectrum`.
"""
function detect(emitted::Dict{<:XRay, <:Real}, det::EDSDetector, response::Matrix{Float64}; noise=false)
    data = zeros(Float64, channelcount(det))
    for (xr, i) in emitted
        ch = channel(energy(xr), det)
        if checkbounds(Bool, data, ch)
            @inbounds data[ch] += i
        end
    end
    meas = response * data # As a Float64[]
    if noise 
        meas = map(c->rand(Poisson(c)), meas)  # Returns an Int[]
    end
    return Spectrum(det.scale, meas, Dict{Symbol,Any}(:Detector=>det))
end

""""
  simulate(comp::Material, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64, det::Detector, resp::Matrix{Float64}; noise=false, vargs...)

Compute a simulated X-ray spectrum for the specified composition material.

Arguments:

  * comp: The Material
  * dose: The electron dose in nA⋅s 
  * e0:   The beam energy in eV
  * θtoa: The take-off angle in radians
  * Ω:    The detector solid angle in steradians
  * det:  The detector model
  * resp: The detector response

Returns a `Spectrum` struct.
"""
function simulate(comp::Material, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64, det::Detector, resp::Matrix{<:AbstractFloat}; noise=false, vargs...)
    ei = emitted_intensities(comp, dose, e0, θtoa, Ω; vargs...)
    sp = detect(ei, det, resp, noise=noise)
    sp[:TakeOffAngle] = θtoa
    sp[:SolidAngle] = Ω
    sp[:ProbeCurrent] = 1.0
    sp[:LiveTime] = dose
    sp[:Composition] = comp
    sp[:BeamEnergy] = e0
    sp[:Name] = "Simulated $(name(comp))"
    return sp
end
