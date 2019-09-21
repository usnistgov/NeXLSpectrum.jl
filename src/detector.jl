using Polynomials


"""
    EnergyScale

An EnergyScale is a way of representing the energy axis associated with X-ray
data. The scale may be linear, polynomial or ??? to handle the various different
non-linearities that happen with EDS detectors plus we can also handle WDS
wavescans.

Implements:

    channel(eV::AbstractFloat, sc::EnergyScale)::Int
    energy(ch::Int, sc::EnergyScale)::Float64
"""
abstract type EnergyScale end

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
        @assert(width>0.0,"Channel width must be greater than 0 in LinearEnergyScale.")
        return new(off,width)
    end
end

Base.show(io::IO, es::LinearEnergyScale) =
    print(io, "E[ch] = ",es.offset," + ",es.width,"⋅ch")

"""
    channel(eV::AbstractFloat, sc::LinearEnergyScale)

    Returns the integer index of the channel for the specified energy X-ray (in
eV).
"""
channel(eV::AbstractFloat, sc::LinearEnergyScale)::Int =
    1 + floor(Int, (eV - sc.offset)/sc.width)


"""
    energy(ch::Integer, sc::LinearEnergyScale)

Returns the energy (in eV) for the low energy side of the bin representing the
ch-th channel.

Example:

    les = LinearEnergyScale(3.0, 10.1)
    energy(101,lsc) == 10.1*101 + 3.0
    energy(101,lsc) - energy(100,lsc) == 10.1
"""
energy(ch::Int, sc::LinearEnergyScale)::Float64 =
    (ch-1)*sc.width+sc.offset


"""
    PolyEnergyScale

An energy scale based on a polynomial function of the channel index.
"""
struct PolyEnergyScale <: EnergyScale
    poly::Poly
    PolyEnergyScale(vec::Vector) = new(Poly(vec,:ch))
end

function Base.show(io::IO, pes::PolyEnergyScale)
    print(io, "E[ch] = ")
    printpoly(io, pes.poly)
end

"""
    channel(eV::AbstractFloat, sc::PolyEnergyScale)

    Returns the integer index of the channel for the specified energy X-ray (in
eV).
"""
function channel(eV::AbstractFloat, sc::PolyEnergyScale)::Int
    cp=Vector(coeffs(sc.poly))
    cp[1]-=eV
    rts = roots(Poly(cp))
    best = 100000
    for rt in rts
        if !(rt isa Complex)
            best = (rt>=-1000.0) && (rt<best) ? 1+floor(rt) : best
        end
    end
    return best
end


"""
    energy(ch::Integer, sc::PolyEnergyScale)

Returns the energy (in eV) for the low energy side of the bin representing the
ch-th channel.

Example:

    pes = PolyEnergyScale([ 3.0, 10.1, 0.001])
    energy(101,pes) ==  3.0 + 10.0*101 + 0.001*101^2
"""
energy(ch::Integer, sc::PolyEnergyScale)::Float64 =
    sc.poly(convert(Float64,ch-1))


"""
    energyscale(es::EnergyScale, channels::UnitRange{Int})

Computes the energy associated with a range of channel indexes and returns
it as an Array.
"""
energyscale(es::EnergyScale, channels::UnitRange{Int}) =
    collect(energy(ch, es) for ch in channels)

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
# Implements linewidth(eV::Float64, res::Resolution)::Float

struct MnKaResolution <: Resolution
    fwhmatmnka::Float64
end

""""
    linewidth(eV::Float64, fwhm::Float64, fwhmenergy::Float64 = 5898.7)

Chuck Fiori's simple function relating the FWHM at eV to the FWHM at another energy.
"""
resolution(eV::Float64, res::MnKaResolution) =
    sqrt((2.45 * (eV - 5898.7)) + (res.fwhmatmnka * res.fwhmatmnka))
"""
    gaussianwidth(fwhm::Float64)

Converts full-width half-max to Gaussian width.
"""
gaussianwidth(fwhm::Float64) =
    fwhm / 2.354820045030949382023138652918

""""
    profile(energy::Float64, xrayE::Float64, res::MnKaResolution)
Calculates a Gaussian profile for an X-ray of xrayE (eV) for a detector
with the specified resolution.
"""
profile(energy::Float64, xrayE::Float64, res::MnKaResolution) =
    exp(-0.5*((energy-xrayE)/gaussianwidth(resolution(xrayE, res)))^2)

"""
    extent(xrayE::Float64, res::MnKaPlusICC, ampl::Float64)
Calculates the extent of the peak interval for an x-ray of the specified
energy.
"""
function extent(xrayE::Float64, res::MnKaResolution, ampl::Float64)
    w=gaussianwidth(resolution(xrayE, res))*sqrt(-2.0*log(ampl))
    return (xrayE-w, xrayE+w)
end

struct MnKaPlusICC <: Resolution
    fwhmatmnka::Float64
    icc::Float64
    iccmax::Float64
    MnKaPlusICC(fwhm::Float64, icc::Float64=70.0, iccmax::Float64=1500.0) = new(fwhm, icc, iccmax)
end

""""
    linewidth(eV::Float64, fwhm::Float64, fwhmenergy::Float64 = 5898.7)

Chuck Fiori's simple function relating the FWHM at eV to the FWHM at another energy plus
a term to account for incomplete charge collection.
"""
resolution(eV::Float64, res::MnKaPlusICC) =
    sqrt((2.45 * (eV - 5898.7)) + (res.fwhmatmnka * res.fwhmatmnka)) + max(0.0, res.icc*(1.0-max(0.0,eV)/res.iccmax))

""""
    profile(energy::Float64, xrayE::Float64, res::MnKaResolution)

Calculates a Gaussian profile for an X-ray of xrayE (eV) for a detector
with the specified resolution.
"""
profile(energy::Float64, xrayE::Float64, res::MnKaPlusICC) =
    exp(-0.5*((energy-xrayE)/gaussianwidth(resolution(xrayE, res)))^2)


"""
    extent(xrayE::Float64, res::MnKaPlusICC, ampl::Float64)

Calculates the extent of the peak interval for an x-ray of the specified
energy.  This includes extra at the low-energy side to account for ICC.
"""
function extent(xrayE::Float64, res::MnKaPlusICC, ampl::Float64)
    w=gaussianwidth(resolution(xrayE, res))*sqrt(-2.0*log(ampl))
    icc=gaussianwidth(res.icc*(1.0-max(0.0,xrayE)/res.iccmax))
    return (xrayE-(w+max(0.0,icc)), xrayE+w)
end

"""
    extent(cxrs::AbstractArray{CharXRay}, res::Resolution, ampl::Float64)::Tuple{2,Float64}

Computes the energy range encompassed by the specified set of x-ray transitions
down to an intensity of ampl.  Relative line weights are taken into account.
"""
function extent(cxrs::AbstractArray{CharXRay,1}, res::Resolution, ampl::Float64)
    tmp = map(cxr -> extent(energy(cxr) ,res, min(0.999,ampl/weight(cxr))), cxrs)
    low, high = +100000,-100000
    for (l, h) in tmp
        low = min(l, low)
        high = max(h, high)
    end
    return (low, high)
end


"""
    Detector

An abstract type defining the characteristics of an X-ray detector.

Implements:

    channelcount(det::Detector)::Int
    scale(det::Detector)::EnergyScale
    resolution(eV::Float64, det::Detector)::Float64
    energy(ch::Int, det::Detector)::Float64
    channel(eV::Float64, det::Detector)::Int
    profile(energy::Float64, xrayE::Float64, det::Detector)
"""
abstract type Detector end


"""
    SimpleEDS

An implementation of the Detector abstract structure for a basic EDS detector.
It includes models of the number of channels, EnergyScale and Resolution.
"""
struct SimpleEDS <: Detector
    channelcount::Int
    scale::EnergyScale
    resolution::Resolution
end

"""
    channelcount(det::SimpleEDS)

Number of detector channels.
"""
channelcount(det::SimpleEDS) =
    det.channelcount

"""
    scale(det::SimpleEDS)

EnergyScale associated with this detector.
"""
scale(det::SimpleEDS)::EnergyScale =
    det.scale

"""
    resolution(eV::Float64, det::SimpleEDS)

Resolution of the detector in eV at the specified energy.
"""
resolution(eV::Float64, det::SimpleEDS) =
    resolution(eV, det.resolution)


"""
    energy(ch::Int, det::SimpleEDS)

Energy of the low-energy side of the ch-th detector bin.
"""
energy(ch::Int, det::SimpleEDS) =
    energy(ch, det.scale)

"""
    channel(eV::Float64, det::SimpleEDS)

The channel index in which the specified energy X-ray belongs.
"""
channel(eV::Float64, det::SimpleEDS) =
    channel(eV, det.scale)

""""
    profile(energy::Float64, xrayE::Float64, det::Detector)

Calculates the profile for an X-ray of xrayE (eV) for a detector
with the specified resolution.
"""
profile(energy::Float64, xrayE::Float64, det::Detector) =
    profile(energy, xrayE, det.resolution)

"""
    extent(cxrs::AbstractArray{CharXRay}, det::Detector, ampl::Float64)::UnitRange

Computes the channel range encompassed by the specified set of x-ray transitions
down to an intensity of ampl.  Relative line weights are taken into account.
"""
extent(cxrs::AbstractArray{CharXRay}, det::Detector, ampl::Float64) =
    UnitRange(map(ee->channel(ee,det), extent(cxrs, det.resolution, ampl))...)

"""
    basicEDS(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64)

Construct simple model of an EDS detector.
"""
basicEDS(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64) =
    SimpleEDS(chCount, LinearEnergyScale(offset,width), MnKaResolution(fwhmatmnka))

"""
    basicEDSwICC(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64)

Construct simple model of an EDS detector with incomplete charge collection at
low X-ray energies.
"""
basicEDSwICC(chCount::Integer, width::Float64, offset::Float64, fwhmatmnka::Float64) =
    SimpleEDS(chCount, LinearEnergyScale(offset,width), MnKaPlusICC(fwhmatmnka, 70.0, 1200.0))
