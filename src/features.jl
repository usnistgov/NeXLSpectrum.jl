
# Features are peaks the may occur in a spectrum. They include characteristic X-rays,
# Compton shifted characteristic X-rays, and escape peaks.

struct EscapeArtifact
    xray::CharXRay
    escape::CharXRay
    function EscapeArtifact(xray::CharXRay, escape::CharXRay = n"Si K-L3")
        @assert energy(xray) > energy(escape) "The energy($xray) must be larger than energy($escape)"
        return new(xray, escape)
    end
end

NeXLCore.exists(::Type{EscapeArtifact}, cxr::CharXRay, esc::CharXRay = n"Si K-L3") =
    energy(cxr) - energy(esc) > 100.0

NeXLCore.element(esc::EscapeArtifact) = element(esc.xray)
NeXLCore.brightest(escs::AbstractVector{EscapeArtifact}) =
    brightest(map(ea -> ea.xray, escs))

Base.show(io::IO, esc::EscapeArtifact) = print(io, name(esc))
Base.show(io::IO, escs::AbstractVector{EscapeArtifact}) = print(io, name(escs))
NeXLCore.name(esc::EscapeArtifact) = "Esc[$(esc.xray)]"
NeXLCore.name(escs::AbstractVector{EscapeArtifact}) =
    "Esc[$(join(map(esc->"$(esc.xray)", escs),","))]"


struct ComptonArtifact
    xray::CharXRay
    angle::Float64
    function ComptonArtifact(prim::CharXRay, θ::Float64)
        @assert (θ >= 0.0) && (θ <= π) "The Compton scatter angle must be between 0.0 and π"
        return new(prim, θ)
    end
end

Base.show(io::IO, ca::ComptonArtifact) = print(io, name(ca))
NeXLCore.name(ca::ComptonArtifact) = "Compton[$(ca.xray), $(rad2deg(ca.angle))°]"
NeXLCore.name(cas::AbstractVector{ComptonArtifact}) =
    "Compton[$(name(map(ca->ca.xray, cas)))]"
NeXLCore.element(ca::ComptonArtifact) = element(xray)

NeXLCore.energy(ca::ComptonArtifact) =
    energy(ca.xray) * NeXLCore.comptonShift(ca.angle, energy(ca.xray))
NeXLCore.weight(esc::ComptonArtifact) = weight(NormalizeToUnity, esc.xray)

NeXLCore.energy(esc::EscapeArtifact) = energy(esc.xray) - energy(esc.escape)

"""
    weight(esc::EscapeArtifact, factor=0.01)

The weight of an EscapeArtifact which is factor * weight(esc.xray).
"""
NeXLCore.weight(esc::EscapeArtifact, factor = 0.01) = factor * weight(NormalizeToUnity, esc.xray)


"""
    SpectrumFeature

A union representing the different type of peak-like features (helpful and harmful) that can appear in a spectrum.
"""
SpectrumFeature = Union{CharXRay,EscapeArtifact,ComptonArtifact}

Base.show(io::IO, sfs::AbstractVector{<:SpectrumFeature}) =
    map(sf -> println(io, name(sf)), sfs)

Base.convert(
    ::AbstractVector{SpectrumFeature},
    x::AbstractVector{CharXRay},
)::AbstractVector{SpectrumFeature} = x

charFeature(
    elm::Element,
    tr::Transition...;
    minweight = 1.0e-3,
    maxE = 1.0e6,
)::AbstractVector{SpectrumFeature} = characteristic(elm, tr, minweight, maxE)

escapeFeature(
    elm::Element,
    trs::Transition...;
    minweight = 0.1,
    maxE = 1.0e6,
    escape = n"Si K-L3",
)::AbstractVector{SpectrumFeature} =
    map(tr -> EscapeArtifact(tr, escape), characteristic(elm, trs, minweight, maxE))

comptonFeature(
    elm::Element,
    trs::Transition...;
    minweight = 1.0e-3,
)::AbstractVector{SpectrumFeature} =
    map(tr -> ComptonArtifact(tr, escape), characteristic(elm, trs, minweight, maxE))
