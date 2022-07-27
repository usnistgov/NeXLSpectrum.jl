"""
    χ²(s1::Spectrum{T}, s2::Spectrum{U}, chs)::T where {T<:Real, U <: Real}
    χ²(specs::AbstractArray{Spectrum{T}}, chs)::Matrix{T}

Computes the dose corrected reduced χ² metric between `s1` and `s2` over the channels in `chs`.

The second form computes a matrix of χ² comparing each spectrum in the array to the others.
"""
function χ²(s1::Spectrum{T}, s2::Spectrum{U}, chs)::Float64 where {T<:Real,U<:Real}
    k1, k2 = 1.0 / dose(s1), 1.0 / dose(s2)
    return sum(chs) do ch
        s1c, s2c = get(s1, ch, 0.0), get(s2, ch, 0.0)
        (k1 * s1c - k2 * s2c)^2 / (k1 * k1 * max(one(T), s1c) + k2 * k2 * max(one(U), s2c))
    end
end
function χ²(specs::AbstractVector{<:Spectrum}, chs)::Matrix{Float64}
    χ2s = zeros(Float64, (length(specs), length(specs)))
    for i in eachindex(specs), j in i+1:length(specs)
        χ2s[i, j] = χ²(specs[i], specs[j], chs)
        χ2s[j, i] = χ2s[i, j]
    end
    return χ2s
end

"""
    similarity(s1::Spectrum{T}, s2::Spectrum{T}, chs)::Float64 where {T<:Real}
    similarity(specs::AbstractArray{Spectrum{T}}, chs)::Vector{Float64}
    similarity(specs::AbstractArray{Spectrum}, minE::Float64=100.0)::Vector{Float64}
    similarity(specs::AbstractArray{<:Spectrum}, det::Detector, elm::Element)::Vector{Float64}
    similarity(specs::AbstractArray{Spectrum{T}}, det::Detector, mat::Material)::Vector{Float64}

Returns a vector of similarity metrics which measure how similar the i-th `Spectrum` is to the other spectra.
The mean reduced χ² statistic metric is such that if `s1` and `s2` differ by only count statistics 
then the metric will be approximately unity.  If `s1` and `s2` vary due to probe current drift, sample 
inhomogeneity, surface roughness or other non-count statistics related reasons then the metric will be 
larger than one.

The first version covers all the channels between minE and the nominal beam energy. The third and fourth versions
considers those channels representing peaks in a spectrum from the `Material` or `Element` on the `Detector`.
"""
similarity(s1::Spectrum, s2::Spectrum, chs)::Float64 = χ²(s1, s2, chs) / length(chs)
function similarity(
    specs::AbstractArray{<:Spectrum},
    chs
)::Vector{Float64}
    return [similarity(spec, sum(filter(s -> !(s === spec), specs)), chs) for spec in specs]
end

function similarity(
    specs::AbstractArray{<:Spectrum},
    minE::Float64=100.0,
)::Vector{Float64}
    e0 = maximum(spec[:BeamEnergy] for spec in specs)
    chs =
        minimum(
            channel(minE, spec) for spec in specs
        ):maximum(channel(e0, spec) for spec in specs)
    return similarity(specs, chs)
end
function similarity(
    specs::AbstractArray{<:Spectrum},
    det::Detector,
    elm::Element,
)::Vector{Float64}
    e0 = maximum(spec[:BeamEnergy] for spec in specs)
    rois = extents(characteristic(elm, alltransitions, 0.01, e0), det, 0.001)
    chs = mapreduce(collect, append!, rois)
    return similarity(specs, chs)
end
function similarity(
    specs::AbstractArray{<:Spectrum},
    det::Detector,
    mat::Material,
)::Vector{Float64}
    function mrg(inp::Vector{UnitRange{Int}})::Vector{UnitRange{Int}}
        simp = sort(inp)
        st, res = simp[1], UnitRange{Int}[]
        for r in simp[2:end]
            if isempty(intersect(st, r))
                push!(res, st)
                st = r
            else
                st = first(st):max(last(r), last(st))
            end
        end
        push!(res, st)
        return res
    end
    e0 = maximum(spec[:BeamEnergy] for spec in specs)
    # Figure out the contiguous ROIs and the channels in the ROIs 
    rois = mrg(
        mapreduce(
            elm -> extents(characteristic(elm, alltransitions, 0.01, e0), det, 0.001),
            append!,
            keys(mat),
        ),
    )
    chs = mapreduce(collect, append!, rois)
    return similarity(specs, chs)
end


"""
    findsimilar(
        specs::AbstractArray{Spectrum{T}}; 
        atol = nothing, 
        rtol=1.5, 
        minspecs=3
    )::Vector{Spectrum{T}}
    findsimilar(
        specs::AbstractArray{Spectrum{T}},
        det::Detector,
        elm::Element; 
        atol = nothing, 
        rtol=1.5, 
        minspecs = 3
    )::Vector{Spectrum{T}}


Filters a collection of spectra for the ones most similar to the average by
removing the least similar spectrum sequentially until all the remaining spectra are either:
 
  * less than atol (if atol != nothing)
  * less than rtol * median(others) (if rtol != nothing)

when applying the 'similarity(...)` function to the spectrum and the sum of the other spectra.

This is useful for finding which of a set of replicate spectra are sufficiently similar 
to each other.
"""
function findsimilar(
    specs::AbstractArray{Spectrum{T}};
    atol=nothing,
    rtol=1.5,
    minspecs=3
)::Vector{Spectrum{T}} where {T<:Real}
    function keep(fs, meds)
        return (isnothing(atol) || (fs < atol)) && #
                (isnothing(rtol) || (fs < rtol * meds))
    end
    if length(specs) >= minspecs
        σs = similarity(specs)
        (fmσ, fmi) = findmax(σs)
        # Now perform this recursively until all are within tol or we hit minspecs
        if !keep(fmσ, median(filter(σ -> σ ≠ fmσ, σs)))
            rem = filter(s -> !(s === specs[fmi]), specs)
            return findsimilar(rem, atol=atol, rtol=rtol, minspecs=minspecs)
        else
            return specs
        end
    end
    error("There are not $minspecs spectra which are sufficiently similar to the mean spectrum.")
end
function findsimilar(
    specs::AbstractArray{Spectrum{T}},
    det::Detector,
    elm::Element;
    atol=nothing,
    rtol=1.5,
    minspecs=3
)::Vector{Spectrum{T}} where {T<:Real}
    function keep(fs, meds)
        return (isnothing(atol) || (fs < atol)) && #
                (isnothing(rtol) || (fs < rtol * meds))
    end
    if length(specs) >= minspecs
        σs = NeXLSpectrum.similarity(specs, det, elm)
        (fmσ, fmi) = findmax(σs)
        # Now perform this recursively until all are within tol or we hit minspecs
        if !keep(fmσ, median(filter(σ -> σ ≠ fmσ, σs)))
            rem = filter(s -> !(s === specs[fmi]), specs)
            return findsimilar(rem, det, elm, atol=atol, rtol=rtol, minspecs=minspecs)
        else
            return specs
        end
    end
    error("There are not $minspecs spectra which are sufficiently similar to the mean spectrum.")
end


"""
    sigma(spec::Spectrum, specs::AbstractArray{<:Spectrum}, chs::AbstractRange{<:Integer})::Vector{Float64}

Computes on a channel-by-channel basis how much `spec` spectrum deviates from the mean of the
other spectra in `specs`.  The result is expressed in terms of the standard deviation expected
from count statistics alone.   Assuming `spec` varies only by count statistics we expect
the result values have a mean 0.0 and a standard deviation of 1.0. 
"""
function sigma(spec::Spectrum, specs::AbstractArray{<:Spectrum}, chs::AbstractRange{<:Integer})::Vector{Float64}
    function doseaverage(specs, chs)
        t = [uv(spec, chs) / dose(spec) for spec in specs]
        return map(j -> mean(collect(t[i][j] for i in eachindex(t))), eachindex(t[1]))
    end
    function delta(spec, specs, chs)
        minus(uv1, uv2) = uv(value(uv1) - value(uv2), sqrt(variance(uv1) + variance(uv2)))
        return minus.(uv(spec, chs), dose(spec) * doseaverage(filter(s -> s != spec, specs), chs))
    end
    return map(v -> value(v) / σ(v), delta(spec, specs, chs))
end

"""
    dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{T} where { T <: Real }
    dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{Float64} where { T <: Integer }

Compute a spectrum which is `spectrum` rescaled to a live time times probe current equal to `dose`.
Useful for setting spectra on an equivalent acquisition duration scale.
"""
function dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{T} where {T<:AbstractFloat}
    res = copy(spectrum)
    scale = dose / NeXLSpectrum.dose(res)
    res.counts .*= scale
    res[:LiveTime] *= scale
    res[:Name] = "N[$(spectrum[:Name]), $dose nA⋅s]"
    return res
end
function dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{Float64} where {T<:Integer}
    scale = dose / NeXLSpectrum.dose(spectrum)
    newProps = copy(spectrum.properties)
    newProps[:LiveTime] *= scale
    newProps[:Name] = "N[$(spectrum[:Name]), $dose nA⋅s]"
    return Spectrum(spectrum.energy, Float64.(spectrum.counts) * scale, newProps)
end


"""
    shannon_entropy(spec::Spectrum)

Computes a measure of the information content in a spectrum.  As there become more and more
distinct values in a spectrum, this value approaches log2(nchannels(spec)).  This number
reflects the number of bits necessary to encode the spectrum data with maximum efficiency.

This is inspired by John Colby's FLAME software which did something similar.  Although, to
be honest, I don't know how his algorithm was implemented.
"""
function shannon_entropy(spec::Spectrum)
    d = Dict{Int,Float64}()
    rr(c::AbstractFloat) = round(Int, c, RoundNearestTiesUp)
    rr(c::Integer) = Int(c)
    foreach(rr.(counts(spec))) do ci
        d[ci] = get(d, ci, 0.0) + 1.0 / length(spec)
    end
    return -sum(cx -> cx * log2(cx), values(d))
end