# Schamber's Vector quant algorithm - It foregoes the a weighted least squares fit for a simple
# homoskedastic fit.  In doing this it can precompute a fitting matrix which requires nothing
# more than a single matrix multiplication to perform the fit.  This makes this mechanism
# extremely quick. This makes it ideal for processing in real-time or HyperSpectrum objects.

struct _VQRefData{T <: AbstractFloat} 
    label::ReferenceLabel
    roi::UnitRange{Int}
    charonly::Vector{T}
    sumchar::T
    scale::T
end

function Base.show(io::IO, vqr::_VQRefData)
    print(io,"$(vqr.label)[$(vqe.roi)]")
end


struct VectorQuant{T <: AbstractFloat} 
    # Vector(label[1], roi[2], charonly[3], sum(charonly)[4], scale[5])
    references::Vector{_VQRefData{T}}
    vectors::Matrix{T}

    """
        VectorQuant(frefs::Vector{FilteredReference}, filt::TopHatFilter)
    
    Constructs a structure used to perform accelerated filtered spectrum fits based on the specified
    collection of `FilteredReference`(s), and a `TopHatFilter`.
    """
    function VectorQuant(frefs::Vector{<:FilteredReference}, filt::TopHatFilter{T}) where { T <: AbstractFloat}
        refs = map(frefs) do fref 
            _VQRefData(fref.label, fref.roi, fref.charonly, sum(fref.charonly), fref.scale)
        end
        x = zeros(T, (length(filt.filters), length(frefs)))
        for (c, fref) in enumerate(frefs)
            x[fref.ffroi, c] = fref.filtered
        end
        # ((ch × ne)T * (ch × ne))^(-1) * (ch × ne) * (ch × ch) => (ne × ch)
        xTxIxf =
            pinv(transpose(x) * x) *
            transpose(x) *
            filterdata(filt, 1:length(filt.filters))
        return new{T}(refs, xTxIxf)
    end

    VectorQuant(ffrs::FilterFitPacket{T}) where { T <: AbstractFloat } = #
        VectorQuant(ffrs.references, ffrs.filter)
end

NeXLCore.minproperties(::VectorQuant) = (:BeamEnergy, :TakeOffAngle, :)

Base.show(io::IO, vq::VectorQuant) = print(
    io,
    "VectorQuant[\n" * join(map(r -> "\t$r\n", vq.references), ",\n") * "\n]",
)

function fit_spectrum(
    vq::VectorQuant{T},
    spec::Spectrum,
    zero = x -> max(0.0, x),
) where { T <: AbstractFloat }
    raw = counts(spec, T)
    krs = zero.(vq.vectors * raw)
    spsc = dose(spec)
    residual = copy(raw)
    for (i, vqr) in enumerate(vq.references)
        residual[vqr.roi] -= krs[i] * vqr.charonly
    end
    peakback = Dict{ReferenceLabel,NTuple{3,Float64}}()
    dkrs = zeros(Float64, length(vq.references))
    for (i, vqr) in enumerate(vq.references)
        ii, bb = krs[i] * vqr.sumchar, sum(residual[vqr.roi])
        peakback[vqr.label] = (ii, bb, bb / spsc)
        dkrs[i] = sqrt(max(0.0, ii + bb)) / vqr.sumchar
    end
    kratios = uvs(
        map(ref -> ref[1], vq.references), #
        map(i -> krs[i] / (vq.references[i].scale * spsc), eachindex(krs)), #
        map(i -> (dkrs[i] / (vq.references[i].scale * spsc))^2, eachindex(krs)),
    )
    return FilterFitResult{T}(
        UnknownLabel(spec),
        kratios,
        1:length(raw),
        raw,
        residual,
        peakback
    )
end

function fit_spectrum(
    vq::VectorQuant,
    hs::HyperSpectrum,
    zero = x -> max(0.0, x),
)::Array{KRatios}
    krs = zeros(Float32, length(vq.references), size(hs)...)
    vecs = vq.vectors[:, 1:depth(hs)]
    scales = [ ref.scale for ref in vq.references ]
    data = counts(hs)
    # @threads seems to slow this (maybe cache misses??)
    for ci in CartesianIndices(hs)
        @inbounds @avx krs[:, ci] = (vecs * data[:, ci]) ./ (dose(hs,ci) * scales)
    end
    # ensure positive...
    map!(zero, krs, krs)
    res = KRatios[]
    for i in filter(ii -> vq.references[ii].label isa CharXRayLabel, eachindex(vq.references))
        lbl = vq.references[i].label
        rprops = properties(spectrum(lbl))
        push!(
            res,
            KRatios(xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i, :, :]),
        )
    end
    return res
end
