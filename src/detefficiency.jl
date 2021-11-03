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
    thickness = 0.0370,
    deadlayer = 30.0e-7,
    entrance = Film(pure(n"Al"), 10.0e-7),
) = DetectorEfficiency(
    "SDD",
    window,
    [entrance, Film(pure(n"Si"), deadlayer)],
    Film(pure(n"Si"), thickness),
)

SiLiEfficiency(
    window::AbstractWindow;
    thickness = 0.250,
    deadlayer = 30.0e-7,
    entrance = Film(pure(n"Al"), 10.0e-7),
) = DetectorEfficiency(
    "Si(Li)",
    window,
    [entrance, Film(pure(n"Si"), deadlayer)],
    Film(pure(n"Si"), thickness),
)

efficiency(aa::DetectorEfficiency, energy::Float64, angle::Float64 = π / 2) =
    transmission(aa.window, energy, angle) *
    (1.0 - transmission(aa.active, energy, angle)) *
    mapreduce(lyr -> transmission(lyr, energy, angle), *, aa.surface)

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
function detectorresponse(
    det::EDSDetector,
    eff::DetectorEfficiency,
    incidence::Float64 = π / 2,
)
    res = zeros(Float64, (channelcount(det), channelcount(det)))
    # An x-ray with energy in ch will be dispersed among a range of channels about ch
    for ch in max(channel(10.0, det), lld(det)):channelcount(det) # X-ray energy by channel
        el, eh = energy(ch, det), energy(ch + 1, det)  # X-ray energy
        fwhm = resolution(0.5 * (eh + el), det)
        roc =
            max(
                lld(det),
                channel(el - 3.0 * fwhm, det),
            ):min(channelcount(det), channel(el + 3.0 * fwhm, det))
        tmp = map(ch2 -> profile(ch2, 0.5 * (eh + el), det), roc)
        effic = 0.5 * (efficiency(eff, eh, incidence) + efficiency(eff, el, incidence))
        res[ch, roc] = (effic/sum(tmp))*tmp # Ensure exact normalization
    end
    return res
end
