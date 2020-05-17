# Schamber's Vector quant algorithm - It foregoes the a weighted least squares fit for a simple
# homosketastic fit.  In doing this it can precompute a fitting matrix which requires nothing
# more than a single matrix multiplication to perform the fit.  This makes this mechanism
# extremely quick. This makes it ideal for processing in real-time or HyperSpectrum objects.
using NeXLSpectrum
using LinearAlgebra

defaspure(c::Material,cxr::CharXRay,e0::Float64,toa::Float64) =  c[element(cxr)]
unityaspure(c::Material,cxr::CharXRay,e0::Float64,toa::Float64) = 1.0

struct VectorQuant
    # Vector(label[1], roi[2], charonly[3], sum(charonly)[4], scale[5])
    references::Vector{Tuple{ReferenceLabel, UnitRange, Vector{Float64}, Float64, Float64}}
    vectors::Matrix{Float64}
"""
    VectorQuant(frefs::Vector{FilteredReference}, filt::TopHatFilter, aspure=defaspure)

where

    aspure(c::Material,cxr::CharXRay,e0::Float64,toa::Float64)

Constructs a structure used to perform accelerated filtered spectrum fits based on the specified
collection of `FilteredReference`(s), a `TopHatFilter` and a composition correction function `aspure(...)`.
The composition correction function is intended to take into account the composition of the reference
by scaling the intensity relative to a pure element standard.

For example, if CaF2 is used as the reference for Ca, we'd scale the intensity to account for the
fact that Ca represents on 0.51333 and that there is a matrix correction factor (ZAF) of approximately 0.9695.
This implies that we need a scale factor of 0.51333*0.9695.   A suitable function is implemented in
NeXLMatrixCorrection.aspure(...).  The default implementation `defaspure(...)`` just scales the intensity to
account for the mass fraction. The alternative `unityaspure(...)` doesn't scale at all.
"""
    VectorQuant(frefs::Vector{FilteredReference}, filt::TopHatFilter, aspure=defaspure) =
        new(_buildVectorQuant(frefs, filt, aspure)...)
end

Base.show(io::IO, vq::VectorQuant) =
    print(io, "VectorQuant[\n"*join(map(r->"\t"*repr(r[1]),vq.references),",\n")*"\n]")

function _buildVectorQuant(frefs::Vector{FilteredReference}, filt::TopHatFilter, aspure::Function)
    scale(cxrlbl::CharXRayLabel) = #
        aspure(cxrlbl.spec[:Composition], brightest(cxrlbl.xrays), cxrlbl.spec[:BeamEnergy], cxrlbl.spec[:TakeOffAngle])
    scale(reflbl::ReferenceLabel) = 1.0
    refs = [ ( fref.identifier, fref.roi, fref.charonly, sum(fref.charonly), fref.scale/scale(fref.identifier)) for fref in frefs]
    x = zeros(Float64, (length(filt.filters), length(frefs)))
    for (c, fref) in enumerate(frefs)
        x[fref.ffroi, c] = fref.filtered
    end
    # ((ch × ne)T * (ch × ne))^(-1) * (ch × ne) * (ch × ch) => (ne × ch)
    xTxIxf = pinv(transpose(x)*x)*transpose(x)*NeXLSpectrum.filterdata(filt, 1:length(filt.filters))
    return ( refs, xTxIxf)
end

function NeXLSpectrum.fit(vq::VectorQuant, spec::Spectrum)::FilterFitResult
    raw = counts(spec, Float64)
    krs=vq.vectors*raw
    spsc = dose(spec)
    residual = copy(raw)
    for (i, (_, roi, co, _, _)) in enumerate(vq.references)
        residual[roi] -= krs[i] * co
    end
    peakback = Dict{ReferenceLabel,NTuple{2,Float64}}()
    dkrs = zeros(Float64, length(vq.references))
    for (i, (lbl, roi, _, ico, _)) in enumerate(vq.references)
        ii, bb = krs[i]*ico, sum(residual[roi])
        peakback[lbl]= ( ii, bb)
        dkrs[i]=sqrt(ii+bb)/ico
    end
    kratios = uvs(map(ref->ref[1], vq.references), #
        map(i->krs[i]/(vq.references[i][5]*spsc), eachindex(krs)), #
        map(i->(dkrs[i]/(vq.references[i][5]*spsc))^2, eachindex(krs)))
    return FilterFitResult(UnknownLabel(spec), kratios, 1:length(raw), raw, residual, peakback)
end

function NeXLSpectrum.fit(vq::VectorQuant, hs::HyperSpectrum)
    res = zeros(Float64, ( size(hs)..., length(vq.references)))
    for idx in eachindex(hs)
        rawks=vq.vectors*counts(hs[idx], Float64)
        res[idx] = rawks ./ sum(rawks)
    end
    return ( map(r->r[1], vq.references), res )
end
