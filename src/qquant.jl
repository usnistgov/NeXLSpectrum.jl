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
    print(io,"$(vqr.label)[$(vqr.roi)]")
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

    VectorQuant(ffrs::FilterFitPacket{S, T}) where { S <: Detector, T <: AbstractFloat } = #
        VectorQuant(ffrs.references, ffrs.filter)
end

NeXLCore.minproperties(::VectorQuant) = (:BeamEnergy, :TakeOffAngle, :)

Base.show(io::IO, vq::VectorQuant) = print(
    io,
    "VectorQuant[\n" * join(map(r -> "\t$r", vq.references), ",\n") * "\n]",
)
"""
    fit_spectrum(
        hs::Spectrum|HyperSpectrum,
        vq::VectorQuant{S <: Detector, T <: AbstractFloat},
        zero = x -> max(Base.zero(T), x)
    )

Fit the spectrum or hyper-spectrum using the vector-quant algorithm. The function `zero` is
applied to the resultant k-ratios before they are returned.  The default simply sets negative
k-ratios to 0.0.  `zero=identity` would leave the negative k-ratios as such.
"""
function fit_spectrum(
    spec::Spectrum,
    vq::VectorQuant{T},
    zero = x -> max(Base.zero(T), x),
) where { T <: AbstractFloat }
    raw = counts(spec, 1:size(vq.vectors, 2), T, true)
    krs = zero.(vq.vectors * raw)
    spsc = T(dose(spec))
    residual = copy(raw)
    for (i, vqr) in enumerate(vq.references)
        residual[vqr.roi] -= krs[i] * vqr.charonly
    end
    peakback = Dict{ReferenceLabel,NTuple{3,T}}()
    dkrs = zeros(T, length(vq.references))
    for (i, vqr) in enumerate(vq.references)
        ii, bb = krs[i] * vqr.sumchar, sum(residual[vqr.roi])
        peakback[vqr.label] = (ii, bb, bb / spsc)
        dkrs[i] = sqrt(max(Base.zero(T), ii + bb)) / vqr.sumchar
    end
    kratios = uvs(
        map(ref -> ref.label, vq.references), #
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
function fit_spectra(
    hs::HyperSpectrum,
    vq::VectorQuant{T},
    zero = x -> max(Base.zero(T), x),
)::Array{KRatios} where { T <: AbstractFloat }
    krs = zeros(T, length(vq.references), size(hs)...)
    vecs = vq.vectors[:, 1:depth(hs)]
    foreach(i->vecs[i,:]/=vq.references[i].scale, eachindex(vq.references))
    data = T.(counts(hs))  # One large allocation over many smaller????
    # In my testing, ThreadsX produces about a 2.7x speedup for 4 threads on 4 cores
    ThreadsX.foreach(CartesianIndices(hs)) do ci
        krs[:, ci] .= (vecs * view(data, :, ci)) / T(dose(hs, ci))
    end
    # ensure positive...
    map!(zero, krs, krs)
    res = map(filter(ii -> vq.references[ii].label isa CharXRayLabel, eachindex(vq.references))) do i
        lbl = vq.references[i].label
        rprops = properties(spectrum(lbl))
        KRatios(xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i, :, :])
    end
    return res
end

fit_spectrum(
    hs::HyperSpectrum,
    vq::VectorQuant{T},
    zero = x -> max(Base.zero(T), x),
) where { T<: AbstractFloat } = fit_spectra(hs, vq, zero)