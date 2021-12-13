
"""
    EnergyScale

An EnergyScale is a way of representing the energy axis associated with X-ray
data. The scale may be linear, polynomial or ??? to handle the various different
non-linearities that happen with EDS detectors plus we can also handle WDS
wavescans.

Implements:

    channel(::Type{Float64}, eV::AbstractFloat, sc::EnergyScale)::Float64
    energy(ch::Int, sc::EnergyScale)::Float64
"""
abstract type EnergyScale end

channel(eV::AbstractFloat, sc::EnergyScale)::Int =
    1 + floor(Int, channel(Float64, eV, sc))

"""
    LinearEnergyScale

An EnergyScale implementation parameterized by channel width and offset.
"""
struct LinearEnergyScale <: EnergyScale
    offset::Float64
    width::Float64
    """
        LinearEnergyScale(off::Float64, width::Float64)
    
    Construct a LinearEnergyScale from the zero offset (in eV) and channel width (in
    eV/channel).
    """
    function LinearEnergyScale(off::Float64, width::Float64)
        @assert width > 0.0 "Channel width must be greater than 0 in LinearEnergyScale."
        return new(off, width)
    end
end

Base.show(io::IO, es::LinearEnergyScale) =
    print(io, es.offset, " + ", es.width, "⋅ch eV")

Base.hash(les::LinearEnergyScale, h::UInt64) = hash(les.offset, hash(les.width, h))

"""
    channel(::Type{Float64}, eV::AbstractFloat, sc::EnergyScale)::Float64

Returns the fractional-channel index of the specified energy X-ray (in eV).

    channel(eV::AbstractFloat, sc::EnergyScale)::Int
    channel(eV::Float64, spec::Spectrum)::Int
    channel(eV::Float64, det::EDSDetector)::Int

Returns the integer index of the channel for the specified energy X-ray (in eV).
"""
channel(::Type{Float64}, eV::AbstractFloat, sc::LinearEnergyScale) = (eV - sc.offset) / sc.width


"""
    NeXLCore.energy(ch::Integer, sc::EnergyScale)
    NeXLCore.energy(ch::Int, spec::Spectrum)
    NeXLCore.energy(ch::Int, det::EDSDetector)


Returns the energy (in eV) for the low energy side of the bin representing the
ch-th channel.

Example:

    les = LinearEnergyScale(3.0, 10.1)
    energy(101,lsc) == 10.1*101 + 3.0
    energy(101,lsc) - energy(100,lsc) == 10.1
    pes = PolyEnergyScale([ 3.0, 10.1, 0.001])
    energy(101,pes) ==  3.0 + 10.0*101 + 0.001*101^2
"""
NeXLCore.energy(ch::Int, sc::LinearEnergyScale)::Float64 = (ch - 1) * sc.width + sc.offset


"""
    PolyEnergyScale

An energy scale based on a polynomial function of the channel index.
"""
struct PolyEnergyScale <: EnergyScale
    poly::ImmutablePolynomial
    PolyEnergyScale(vec::AbstractVector{<:AbstractFloat}) =
        new(ImmutablePolynomial(vec, :ch))
end

function Base.show(io::IO, pes::PolyEnergyScale)
    printpoly(io, pes.poly)
    print(io, " eV")
end

function channel(::Type{Float64}, eV::AbstractFloat, sc::PolyEnergyScale)::Float64
    cp = Float64[ coeffs(sc.poly)... ]
    cp[1] -= eV
    rts = roots(ImmutablePolynomial(cp))
    best = 100000.0
    for rt in rts
        if abs(imag(rt))<1.0e-4
            rrt = real(rt)
            best = (rrt >= -1000.0) && (rrt < best) ? rrt : best
        end
    end
    return best
end

NeXLCore.energy(ch::Integer, sc::PolyEnergyScale)::Float64 =
    sc.poly(convert(Float64, ch - 1))


"""
    energyscale(es::EnergyScale, channels)
    energyscale(det::EDSDetector)

Computes the energy associated with a range of channel indexes and returns
it as an Array.
"""
energyscale(es::EnergyScale, channels) = map(ch->energy(ch, es), channels)

"""
    Resolution

An abstract type describing the channel dependence of the resolution of an EDS
detector.

Implements:

    resolution(eV::Float64, res::Resolution)::Float # Resolution at specified energy
    profile(energy::Float64, xrayE::Float64, res::Resolution) # Amplitude for a signal at the specified energy at the specified energy
    extent(xrayE::Float64, res::Resolution, ampl::Float64)::Tuple{2,Float} # The range of channels over which the signal exceeds ampl
"""
abstract type Resolution end

"""
    MnKaResolution

Uses Chuck Fiori's simple function relating the FWHM at eV to the FWHM at another energy.
"""
struct MnKaResolution <: Resolution
    fwhmatmnka::Float64
end

Base.show(io::IO, mnka::MnKaResolution) = print(io, "$(mnka.fwhmatmnka) eV @ Mn K-L3")


"""
    resolution_to_fwhm(::Type{MnKaResolution}, res::Float64, eV::Float64)

Given the FWHM at res predict the resolution at Mn Kα.
"""
resolution_to_fwhm(::Type{MnKaResolution}, res::Float64, eV::Float64) = 
    sqrt(2.45 * (eV - enx"Mn K-L3") + res^2)


""""
    resolution(eV::Float64, res::Resolution)
    resolution(eV::Float64, det::EDSDetector)

The FWHM at eV for the `<:Resolution` model.
"""
resolution(eV::Float64, res::MnKaResolution) = resolution_to_fwhm(MnKaResolution, res.fwhmatmnka, eV)

"""
    gaussianwidth(fwhm::Float64)

Converts full-width half-max to Gaussian width.  See also fwhm(...)
"""
gaussianwidth(fwhm::Float64) = fwhm / (2.0*sqrt(2.0*log(2.0)))
"""
    fwhm(gauss::Float64)

    
Converts Gaussian width to full-width half-max.  See also gaussianwidth
"""
fwhm(gauss::Float64) = gauss*(2.0*sqrt(2.0*log(2.0)))

""""
    profile(energy::Float64, xrayE::Float64, res::Resolution)

Calculates a Gaussian profile for an X-ray of xrayE (eV) for a detector
with the specified resolution.  Maintains normalization to a sum of unity.
"""
function profile(energy::Float64, xrayE::Float64, res::MnKaResolution)
    σ = gaussianwidth(resolution(xrayE, res))
    return exp(-0.5 * ((energy - xrayE) / σ)^2) / (σ * sqrt(2.0π))
end

"""
    extent(xrayE::Float64, res::Resolution, ampl::Float64)

Calculates the extent of the peak interval for an x-ray of the specified
energy.
"""
function extent(xrayE::Float64, res::MnKaResolution, ampl::Float64)
    w = gaussianwidth(resolution(xrayE, res)) * sqrt(-2.0 * log(ampl))
    return (xrayE - w, xrayE + w)
end

struct MnKaPlusICC <: Resolution
    fwhmatmnka::Float64
    icc::Float64
    iccmax::Float64
    MnKaPlusICC(fwhm::Float64, icc::Float64 = 70.0, iccmax::Float64 = 1500.0) =
        new(fwhm, icc, iccmax)
end

resolution(eV::Float64, res::MnKaPlusICC) =
    sqrt((2.45 * (eV - 5898.7)) + (res.fwhmatmnka * res.fwhmatmnka)) +
    max(0.0, res.icc * (1.0 - max(0.0, eV) / res.iccmax))

function profile(energy::Float64, xrayE::Float64, res::MnKaPlusICC)
    σ = gaussianwidth(resolution(xrayE, res))
    return exp(-0.5 * ((energy - xrayE) / σ)^2) / (σ * sqrt(2.0π))
end

function extent(xrayE::Float64, res::MnKaPlusICC, ampl::Float64)
    w = gaussianwidth(resolution(xrayE, res)) * sqrt(-2.0 * log(ampl))
    icc = gaussianwidth(res.icc * (1.0 - max(0.0, xrayE) / res.iccmax))
    return (xrayE - (w + max(0.0, icc)), xrayE + w)
end

extent(cxr::CharXRay, res::Resolution, ampl::Float64) =
    extent(energy(cxr), res, min(0.999, ampl / weight(cxr)))

"""
    Detector

An abstract type defining the characteristics of an X-ray detector.

Implements:

    channelcount(det::Detector)::Int
    scale(det::Detector)::EnergyScale
    resolution(eV::Float64, det::Detector)::Float64 # FWHM at eV
    energy(ch::Int, det::Detector)::Float64
    channel(eV::Float64, det::Detector)::Int
    profile(energy::Float64, xrayE::Float64, det::Detector)
    lld(det::Detector)::Int
    isvisible(sf::SpectrumFeature, det::Detector)
"""
abstract type Detector end

"""
    extent(escape::EscapeArtifact, res::Resolution, ampl::Float64)::Tuple{2,Float64}

The extent of an escape artifact is determined by the resolution of the detector at the energy of the escape peak.
"""
extent(esc::EscapeArtifact, res::Resolution, ampl::Float64 = 0.01) =
    extent(energy(esc), res, min(0.999, ampl / weight(esc.xray)))

"""
    extent(escape::ComptonArtifact, res::Resolution, ampl::Float64)::Tuple{2,Float64}

The extent of a Compton artifact is determined by the resolution of the detector at the energy of the Compton peak.
"""
extent(ca::ComptonArtifact, res::Resolution, ampl::Float64 = 1.0e-4) =
    extent(energy(ca), res, min(0.999, ampl / weight(ca.xray)))

"""
    extent(sf::SpectrumFeature, det::Detector, ampl::Float64)::Tuple{Float64, Float64}

Computes the channel range encompassed by the specified set of x-ray transitions
down to an intensity of ampl.  Relative line weights are taken into account.
"""
extent(sf::SpectrumFeature, det::Detector, ampl::Float64)::Tuple{Float64,Float64} =
    extent(sf, det.resolution, ampl)

"""
    EDSDetector

Types extending EDSDetector must have member variables

    channelcount::Int # Number of channels
    scale::EnergyScale # Detector calibration funtion
    resolution::Resolution # Detector lineshape function
    lld::Int # low level discriminator
"""
abstract type EDSDetector <: Detector end

"""
    channelcount(det::Detector)

Number of detector channels.
"""
channelcount(det::EDSDetector) = det.channelcount

Base.eachindex(det::EDSDetector) = Base.OneTo(det.channelcount)

"""
    scale(det::Detector)

EnergyScale associated with this detector.
"""
scale(det::EDSDetector)::EnergyScale = det.scale

energyscale(det::EDSDetector) = energyscale(det.scale, eachindex(det))

"""
    lld(det::EDSDetector)

Low level detection limit in channels.  Channels at or below this value will be zeroed when the lld is applied.
"""
lld(det::EDSDetector) = det.lld

resolution(eV::Float64, det::EDSDetector) = resolution(eV, det.resolution)
resolution(det::Detector) = resolution(enx"Mn K-L3", det)


NeXLCore.energy(ch::Int, det::EDSDetector) = energy(ch, det.scale)

channel(eV::Float64, det::EDSDetector) = channel(eV, det.scale)
channel(sf::SpectrumFeature, det::EDSDetector) = channel(energy(sf), det.scale)
channel(::Type{Float64}, eV::Float64, det::EDSDetector) = channel(Float64, eV, det.scale)
channel(::Type{Float64}, sf::SpectrumFeature, det::EDSDetector) = channel(Float64, energy(sf), det.scale)

channelwidth(ch::Int64, det::EDSDetector) = energy(ch + 1, det) - energy(ch, det)
channelwidth(det::EDSDetector) = (energy(channelcount(det), det) - energy(1, det))/channelcount(det)


""""
    profile(ch::Int, xrayE::Float64, det::EDSDetector)

Calculates the profile for an X-ray of xrayE (eV) for a detector
with the specified resolution.  Performs a crude integration to account for
the channel width.
"""
function profile(ch::Int, xrayE::Float64, det::EDSDetector)
    eh, el = energy(ch + 1, det), energy(ch, det)
    return 0.5 *
           (profile(eh, xrayE, det.resolution) + profile(el, xrayE, det.resolution)) *
           (eh - el)
end
profile(nrg::Float64, xrayE::Float64, det::EDSDetector) =
    profile(channel(nrg, det), xrayE, det.resolution)


"""
    isvisible(cxrs::AbstractVector{<:SpectrumFeature}, det::Detector)

Returns the characteristic x-rays that are visible on the specified detector (ie. Between the LLD and the maximum
channel).
"""
isvisible(sfs::AbstractVector{<:SpectrumFeature}, det::Detector) =
    filter(sf -> isvisible(sf, det), sfs)

"""
    extents(cxrs::AbstractVector{<:SpectrumFeature},det::Detector,ampl::Float64)::Vector{UnitRange{Int}}

Determine the contiguous ranges of channels over which the specified collection of X-rays will be measured on
the specified detector.  The ampl determines the extent of each peak.
"""
function extents( #
    cxrs::AbstractVector{T},
    det::Detector, #
    ampl::Float64,
)::Vector{UnitRange{Int}} where {T<:SpectrumFeature}
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
    inrange(x) = (x.start <= channelcount(det)) && (x.stop > lld(det))
    # Note: extent(cxr,...) takes line weight into account
    es = map(
        cxr -> extent(cxr, det, ampl),
        filter(cxr -> weight(cxr) > ampl, isvisible(cxrs, det)),
    )
    return length(es) > 0 ? filter(
        inrange,
        ascontiguous(map(ee -> channel(ee[1], det):channel(ee[2], det), es)),
    ) : Vector{UnitRange{Int}}[]
end

extents(elm::Element, det::Detector, ampl::Float64)::Vector{UnitRange{Int}} =
    extents(isvisible(characteristic(elm, alltransitions), det), det, ampl)

"""
    function labeledextents(
        cxrs::AbstractVector{T},
        det::Detector,
        ampl::Float64
    )::Vector{Tuple{Vector{T},UnitRange{Int}}} where T <: SpectrumFeature

Creates a vector containing pairs containing a vector of T <: SpectrumFeature and an interval. The interval represents a
contiguous interval over which all the X-rays in the interval are sufficiently close in energy that they will
interfere with each other on the specified detector.
"""
function labeledextents(
    cxrs::AbstractVector{T},
    det::Detector,
    ampl::Float64,
)::Vector{Tuple{Vector{T},UnitRange{Int}}} where {T<:SpectrumFeature}
    fcxrs = filter(cxr -> weight(cxr) > ampl, isvisible(cxrs, det))
    es = map(xr -> extent(xr, det, ampl), fcxrs) # CharXRay -> energy ranges
    le = collect(zip(fcxrs, map(ee -> channel(ee[1], det):channel(ee[2], det), es))) # Energy ranges to channel ranges
    sort!(le, lt = (x1, x2) -> isless(energy(x1[1]), energy(x2[1]))) # sort by x-ray energy
    res = Vector{Tuple{Vector{T},UnitRange{Int}}}()
    if length(le) > 0
        curX, curInt = [le[1][1]], le[1][2]
        for (cxr, interval) in le[2:end]
            if length(intersect(interval, curInt)) > 0 # Add to current extent
                curInt = min(interval.start, curInt.start):max(interval.stop, curInt.stop)
                push!(curX, cxr)
            else # create a new extent
                push!(res, (curX, curInt))
                curX = [cxr]
                curInt = interval
            end
        end
        push!(res, (curX, curInt))
    end
    return res
end

"""
    function escapeextents(
        cxrs::AbstractVector{T},
        det::Detector,
        ampl::Float64
    )::Vector{Tuple{Vector{T},UnitRange{Int}}} where T <: SpectrumFeature

Creates a vector containing pairs containing a vector of T <: SpectrumFeature and an interval. The interval represents a
contiguous interval over which all the X-rays in the interval are sufficiently close in energy that they will
interfere with each other on the specified detector.
"""
function escapeextents(
    cxrs::AbstractVector{T},
    det::Detector,
    ampl::Float64,
    maxE::Float64,
    escape::CharXRay = n"Si K-L3",
    minweight::Float64 = 0.5,
)::Vector{Tuple{Vector{EscapeArtifact},UnitRange{Int}}} where {T<:SpectrumFeature}
    escs = map(
        tr -> EscapeArtifact(tr, escape),
        filter(c -> exists(EscapeArtifact, c, escape), cxrs),
    )
    escs = filter(esc -> isvisible(esc, det), escs)
    es = map(esc -> extent(energy(esc), det.resolution, 0.01), escs) # EscapeFeature -> energy ranges
    le = collect(zip(escs, map(ee -> channel(ee[1], det):channel(ee[2], det), es))) # Energy ranges to channel ranges
    sort!(le, lt = (x1, x2) -> isless(energy(x1[1]), energy(x2[1]))) # sort by x-ray energy
    res = Vector{Tuple{Vector{EscapeArtifact},UnitRange{Int}}}()
    if length(le) > 0
        curX, curInt = [le[1][1]], le[1][2]
        for (cxr, interval) in le[2:end]
            if length(intersect(interval, curInt)) > 0 # Add to current extent
                curInt = min(interval.start, curInt.start):max(interval.stop, curInt.stop)
                push!(curX, cxr)
            else # create a new extent
                push!(res, (curX, curInt))
                curX = [cxr]
                curInt = interval
            end
        end
        push!(res, (curX, curInt))
    end
    return res
end

"""
    simpleEDS(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64, lld::Int = channel(150.0 eV))

Construct simple model of an EDS detector.
"""
function simpleEDS(
    chCount::Integer,
    width::Float64,
    offset::Float64,
    fwhmatmnka::Float64,
    lld::Int = -1,
)
    lld = lld < 1 ? round(Int, (150.0 - offset) / width) : lld
    BasicEDS(chCount, LinearEnergyScale(offset, width), MnKaResolution(fwhmatmnka), lld)
end

"""
    simpleEDSwICC(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64, lld::Int=channel(150.0 eV))

Construct simple model of an EDS detector with incomplete charge collection at
low X-ray energies.
"""
function simpleEDSwICC(
    chCount::Integer,
    width::Float64,
    offset::Float64,
    fwhmatmnka::Float64,
    lld::Int = -1,
)
    lld = lld < 1 ? round(Int, (150.0 - offset) / width) : lld
    BasicEDS(
        chCount,
        LinearEnergyScale(offset, width),
        MnKaPlusICC(fwhmatmnka, 70.0, 1200.0),
        max(1, lld),
    )
end

struct BasicEDS <: EDSDetector
    channelcount::Int
    scale::EnergyScale
    resolution::Resolution
    lld::Int # low level discriminator
    minByShell::Dict{Shell,Element} # Dict( Shell(1)=>n"Be", Shell(2)=>n"Sc", Shell(3)=>n"Ba", Shell(4)=>n"Pu" )
    function BasicEDS(
        channelcount::Int,
        zero::Float64,
        gain::Float64,
        fwhm::Float64,
        lld::Int,
        minByShell::Dict{Shell,Element} = Dict{Shell,Element}(),
    )
        mbs = merge(
            Dict{Shell,Element}(
                KShell => n"Be",
                LShell => n"Sc",
                MShell => n"Ba",
                NShell => n"Pu",
            ),
            minByShell,
        )
        return new(
            channelcount,
            LinearEnergyScale(zero, gain),
            MnKaResolution(fwhm),
            lld,
            mbs,
        )
    end
    function BasicEDS(
        channelcount::Int,
        scale::EnergyScale,
        resolution::Resolution,
        lld::Int,
        minByShell::Dict{Shell,Element} = Dict{Shell,Element}(),
    )
        mbs = merge(
            Dict{Shell,Element}(
                KShell => n"Be",
                LShell => n"Sc",
                MShell => n"Ba",
                NShell => n"Pu",
            ),
            minByShell,
        )
        return new(channelcount, scale, resolution, lld, mbs)
    end
end

function Base.show(io::IO, beds::BasicEDS)
    els = join(symbol.(map(klm -> beds.minByShell[klm], Shell.(1:4))), ",")
    print(
        io,
        "BasicEDS[$(beds.channelcount) chs, $(beds.scale), $(beds.resolution), $(beds.lld) ch LLD, [$els]]",
    )
end

"""
    isvisible(sf::SpectrumFeature, det::Detector)

Is `sf` visible on the specified Detector?
"""
function isvisible(esc::EscapeArtifact, det::BasicEDS)
    ass = AtomicSubShell(det.minByShell[KShell], n"K1")
    return (energy(esc) >= energy(ass)) &&
           (energy(esc) > energy(lld(det), det)) && #
           (energy(esc) < energy(det.channelcount + 1, det))
end

isvisible(cxr::CharXRay, det::BasicEDS) =
    (element(cxr) >= det.minByShell[shell(cxr)]) && #
    (energy(cxr) > energy(lld(det), det)) && #
    (energy(cxr) < energy(det.channelcount + 1, det))
