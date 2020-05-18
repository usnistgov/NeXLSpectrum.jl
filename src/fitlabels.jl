
"""
    FilteredLabel

An abstract type associated with labels of filtered spectrum data objects.  structs that extend FilteredLabel should
have `.spec` members.
"""
abstract type FilteredLabel <: Label end

"""
    spectrum(fl::FilteredLabel)::Spectrum

The spectrum associated with a FilteredLabel-based type.
"""
spectrum(fl::FilteredLabel) = fl.spec

"""
    ReferenceLabel

A label associated with reference spectra.  The label encapsulates the original spectrum and the range of channels
represented by this reference object.  structs that extend ReferenceLabel should have `.roi` and
`.spec` members.
"""
abstract type ReferenceLabel <: FilteredLabel end

Base.isless(rl1::ReferenceLabel, rl2::ReferenceLabel) =
    return isequal(rl1.roi, rl2.roi) ? isless(rl1.spec[:Name], rl2.spec[:Name]) :
           (isequal(rl1.roi.start, rl2.roi.start) ? isless(rl1.roi.stop, rl2.roi.stop) :
            isless(rl1.roi.start, rl2.roi.start))


"""
   channels(rl::ReferenceLabel)::UnitRange{Int}

The range of channels associated with the specified ReferenceLabel.
"""
channels(rl::ReferenceLabel) = rl.roi

Base.show(io::IO, refLab::ReferenceLabel) = print(io::IO, "$(refLab.spec[:Name])[$(refLab.roi)]")

_hashrl(spec,roi,xrays) = xor(xor(hash(spec), hash(roi)), hash(xrays))

Base.isequal(c1::ReferenceLabel, c2::ReferenceLabel) =
    (hash(c1) == hash(c2)) && isequal(c1.roi, c2.roi) && isequal(c1.xrays, c2.xrays) && isequal(c1.spec, c2.spec)

"""
    CharXRayLabel

A ReferenceLabel<:FilteredLabel  that Represents a reference spectrum associated with a set of characteristic x-rays
(CharXRay) objects over a contiguous range of spectrum channels.
"""
struct CharXRayLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
    xrays::Vector{CharXRay}
    hash::UInt
    function CharXRayLabel(spec::Spectrum, roi::UnitRange{Int}, xrays::Vector{CharXRay})
        @assert all(xr->element(xr)==element(xrays[1]), xrays)
        new(spec, roi, xrays, _hashrl(spec,roi,xrays))
    end
end


"""
   xrays(cl::CharXRayLabel)

A list of the X-rays associated with this CharXRayLabel.
"""
xrays(cl::CharXRayLabel) = cl.xrays

NeXLCore.element(cl::CharXRayLabel) = element(cl.xrays[1])

Base.show(io::IO, refLab::CharXRayLabel) = print(io::IO, "$(name(refLab.xrays))")
Base.isequal(rl1::CharXRayLabel, rl2::CharXRayLabel) =
    isequal(rl1.roi, rl2.roi) && isequal(rl1.xrays, rl2.xrays) && isequal(rl1.spec, rl2.spec)

"""
    EscapeLabel

A ReferenceLabel<:FilteredLabel that Represents a reference spectrum associated with an escape peak from a set of
characteristic x-rays (CharXRay) objects over a contiguous range of spectrum channels.
"""
struct EscapeLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
    xrays::Vector{EscapeArtifact}
    hash::UInt

    EscapeLabel(spc::Spectrum, roi::UnitRange{Int}, escs::AbstractVector{EscapeArtifact}) =
        new(spc, roi, convert(Vector{EscapeArtifact}, escs), _hashrl(spc,roi,escs))
end

Base.show(io::IO, escl::EscapeLabel) = print(io, name(escl))
Base.isequal(el1::EscapeLabel, el2::EscapeLabel) =
    isequal(el1.roi, el2.roi) && isequal(el1.xrays, el2.xrays) && isequal(el1.spec, el2.spec)


NeXLCore.name(escl::EscapeLabel) = "Ecs[$(name([esc.xray for esc in escl.xrays]))]"
NeXLCore.element(escl::EscapeLabel) = element(escl.xrays[1])

"""
    UnknownLabel

A FilteredLabel that represents the unknown spectrum.
"""
struct UnknownLabel <: FilteredLabel
    spec::Spectrum
end


Base.show(io::IO, unk::UnknownLabel) = print(io, unk.spec[:Name])
Base.isequal(ul1::UnknownLabel, ul2::UnknownLabel) = isequal(ul1.spec, ul2.spec)
spectrum(unkl::UnknownLabel) = unkl.spec

"""
    HyperSpectrumLabel

A Label that represents a single spectrum with a HyperSpectrum.
"""
struct HyperSpectrumLabel <: Label
    hyperspec::HyperSpectrum
    index::CartesianIndex
    HyperSpectrumLabel(hs::HyperSpectrum, idx::Int...) =
        new(hs, CartesianIndex(idx...))
end

Base.show(io::IO, unk::HyperSpectrumLabel) = print(io, unk.hyperspec[:Name]*"[$(unk.index.I)]")
Base.isequal(ul1::HyperSpectrumLabel, ul2::HyperSpectrumLabel) = (ul1.hyperspec===ul2.hyperspec) && isequal(ul1.index,ul2.index)
spectrum(unkl::HyperSpectrumLabel) = unkl[index]
