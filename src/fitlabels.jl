const SpectrumOrProperties = Union{Spectrum, Dict{Symbol,Any}}

function spectrum(sop::SpectrumOrProperties)::Spectrum 
    @assert sop isa Spectrum "The spectrum is not available."
    sop # Fails fast
end
properties(sop::SpectrumOrProperties)::Dict{Symbol, Any} = sop isa Spectrum ? properties(sop) : sop
hasspectrum(sop::SpectrumOrProperties) = sop isa Spectrum

"""
    ReferenceLabel

A label associated with reference spectra.  The label encapsulates the original spectrum and the range of channels
represented by this reference object.  structs that extend ReferenceLabel should have `.roi`, `.spectrum` and ".hash" 
fields.
"""
abstract type ReferenceLabel <: Label
    #spectrum::SpectrumOrProperties
    #roi::UnitRange{Int}
    #xrays::Vector{CharXRay}
    #hash::UInt
end

"""
   channels(rl::ReferenceLabel)::UnitRange{Int}

The range of channels associated with the specified ReferenceLabel.
"""
channels(rl::ReferenceLabel) = rl.roi
Base.hash(rl::ReferenceLabel, h::UInt) = hash(rl.hash, h)
Base.isequal(c1::ReferenceLabel, c2::ReferenceLabel) =
    (typeof(c1)==typeof(c2)) &&
    (c1.hash == c2.hash) &&
    isequal(c1.roi, c2.roi) &&
    isequal(c1.xrays, c2.xrays) &&
    isequal(c1.spectrum, c2.spectrum)
xrays(cl::ReferenceLabel)::Vector{CharXRay} = cl.xrays
spectrum(cl::ReferenceLabel)::Spectrum = spectrum(cl.spectrum)
properties(cl::ReferenceLabel)::Dict{Symbol,Any} = properties(cl.spectrum)
hasspectrum(cl::ReferenceLabel)::Bool = hasspectrum(cl.spectrum)
composition(cl::ReferenceLabel) = get(properties(cl), :Composition, missing)
NeXLCore.element(cl::ReferenceLabel)::Element = element(cl.xrays[1])
Base.isless(rl1::ReferenceLabel, rl2::ReferenceLabel) =
    return isequal(rl1.roi, rl2.roi) ? isless(rl1.spectrum[:Name], rl2.spectrum[:Name]) :
           (
        isequal(rl1.roi.start, rl2.roi.start) ? isless(rl1.roi.stop, rl2.roi.stop) :
        isless(rl1.roi.start, rl2.roi.start)
    )


"""
    CharXRayLabel

A ReferenceLabel that represents a reference spectrum or reference properties associated with a set of 
characteristic x-rays (CharXRay) objects over a contiguous range of spectrum channels.
"""
struct CharXRayLabel <: ReferenceLabel
    spectrum::SpectrumOrProperties # The spectrum used as the reference for fitting...
    roi::UnitRange{Int}
    xrays::Vector{CharXRay}
    hash::UInt
    function CharXRayLabel(spec::SpectrumOrProperties, roi::UnitRange{Int}, xrays::Vector{CharXRay})
        @assert all(xr -> element(xr) == element(xrays[1]), xrays)
        new(spec, roi, xrays, hash(spec, hash(roi, hash(xrays))))
    end
end

function Base.show(io::IO, cl::CharXRayLabel) 
    comp = composition(cl)
    compname = ismissing(comp) ? "Unspecified" : name(comp)
    print(io,"k[$(name(cl.xrays)), $compname]")
end

"""
    EscapeLabel

A ReferenceLabel<:FilteredLabel that Represents a reference spectrum associated with an escape peak from a set of
characteristic x-rays (CharXRay) objects over a contiguous range of spectrum channels.
"""
struct EscapeLabel <: ReferenceLabel
    spectrum::SpectrumOrProperties
    roi::UnitRange{Int}
    xrays::Vector{EscapeArtifact}
    hash::UInt

    EscapeLabel(spc::SpectrumOrProperties, roi::UnitRange{Int}, escs::AbstractVector{EscapeArtifact}) =
        new(spc, roi, convert(Vector{EscapeArtifact}, escs), hash(spc, hash(roi, hash(escs))))
end

Base.show(io::IO, escl::EscapeLabel) = print(io, name(escl))

NeXLCore.name(escl::EscapeLabel) = "Ecs[$(name([esc.xray for esc in escl.xrays]))]"
NeXLCore.element(escl::EscapeLabel) = element(escl.xrays[1])

"""
    UnknownLabel

A Label that represents the unknown spectrum.
"""
struct UnknownLabel <: Label
    spectrum::Union{HyperSpectrum,Spectrum}
    hash::UInt

    UnknownLabel(spec::Union{HyperSpectrum,Spectrum}) = new(spec, hash(spec))
end

Base.show(io::IO, unk::UnknownLabel) = print(io, unk.spectrum[:Name])
Base.isequal(ul1::UnknownLabel, ul2::UnknownLabel) = (ul1.hash==ul2.hash) && isequal(ul1.spectrum, ul2.spectrum)
spectrum(unkl::UnknownLabel) = unkl.spectrum
properties(unkl::UnknownLabel) = properties(unkl.spectrum)

"""
    HyperSpectrumLabel

A Label that represents a single spectrum with a HyperSpectrum.
"""
struct HyperSpectrumLabel <: Label
    hyperspectrum::HyperSpectrum
    index::CartesianIndex
    HyperSpectrumLabel(hs::HyperSpectrum, idx::Int...) = new(hs, CartesianIndex(idx...))
end

Base.show(io::IO, unk::HyperSpectrumLabel) =
    print(io, unk.hyperspectrum[:Name] * "[$(unk.index.I)]")
Base.isequal(ul1::HyperSpectrumLabel, ul2::HyperSpectrumLabel) =
    (ul1.hyperspectrum === ul2.hyperspec) && isequal(ul1.index, ul2.index)
spectrum(unkl::HyperSpectrumLabel) = unkl.hyperspectrum[unkl.index]
properties(unkl::HyperSpectrumLabel) = properties(unkl.hyperspectrum[unkl.index])
