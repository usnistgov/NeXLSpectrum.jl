# Routines to handle spectra collected simulataneously on multiple detectors

using Formatting

"""
    loadmultispec(path::AbstractString, basefn::AbstractString; indexes=0:3, fnmapper::String = "{1}[{2}].msa")

Load multiple spectra using the `basefn` and `fnmapper` to determine which spectra to load.  The spectra should be
related in the sense that they were all collected simulataneously so they have the same `:RealTime`, `:BeamEnergy` 
and `:LiveTime`.
"""
function loadmultispec(
    path::AbstractString,
    basefn::AbstractString;
    indexes::UnitRange{Int} = 0:3,
    fnmapper::String = "{1}[{2}].msa",
)
    fe = FormatExpr(fnmapper)
    res = [loadspectrum(joinpath(path, format(fe, basefn, i))) for i in indexes]
    checkmultispec(res)
    return res
end

"""
    checkmultispec(specs::AbstractArray{<:Spectrum})

Check that the `:BeamEnergy`, `:RealTime` and `:ProbeCurrent` match for all `specs`.
"""
function checkmultispec(specs::AbstractArray{<:Spectrum})
    s1, so = specs[1], view(specs, 2:length(specs))
    @assert all(s->get(s,:ProbeCurrent,0.0)==get(s1, :ProbeCurrent, 0.0), so) "The spectra must have all been collected at the same probe current."
    @assert all(s->s[:RealTime]==s1[:RealTime], so) "The spectra must have all been collected at the same real time."
    @assert all(s->s[:BeamEnergy]==s1[:BeamEnergy], so) "The spectra must have all been collected at the same beam energy."
    @assert all(s->matches(s, s1), so) "The energy calibration of the spectra must all match (approximately.)"
    @assert all(s->length(s.counts)==length(s1.counts), so) "The channel length of the spectra must all match."
end

"""
   specratio(specs::AbstractArray{<:Spectrum})::Vector{Vector{Float64}}

Computes the channel-by-channel ratio of the counts data over the mean counts data for all spectra for each spectrum.
"""
function specratio(specs::AbstractArray{<:Spectrum})
    ss = map(s -> max(1, s), counts(sum(specs)))
    return length(specs) * [counts(spec) ./ ss for spec in specs]
end
"""
    binnedspecratio(specs::AbstractArray{<:Spectrum}; minE=100.0, deltaE=100.0)::Vector{Vector{Float64}}

Bins the specratio(spec) to reduce the variation.
"""
function binnedspecratio(specs::AbstractArray{<:Spectrum}; minE = 100.0, deltaE = 100.0)
    sr = specratio(specs)
    map(
        i -> map(
            ee -> mean(
                sr[i][channel(
                    ee,
                    specs[i],
                ):channel(min(specs[i][:BeamEnergy], ee + deltaE), specs[i])],
            ),
            minE:deltaE:specs[i][:BeamEnergy]/1.8,
        ),
        eachindex(specs),
    )
end

"""
    multiscore(specs::AbstractArray{<:Spectrum}, e0=specs[1][:BeamEnergy])

Compares each spectrum against the mean spectrum by comparing the 200 eV to 500 eV range with the [e0/2,e0/1.5] range.
If all spectra are equivalent at low energies then all the scores will be close to zero.  A number less than zero
means that the low energy region of the spectrum was depressed relative to the others.  A number more than zero
means that the high energy region of the spectrum was elevated relative to the others.  The multi-score is sensitive
to tilt and obstructions like surface texture which may make one spectrum's low energy be more absorbed than the 
others. 
"""
function multiscore(specs::AbstractArray{<:Spectrum}, e0 = specs[1][:BeamEnergy])
    srs = specratio(specs)
    sp = specs[1]
    rr = channel(200.0, sp):channel(500.0, sp)
    ss = channel(e0 / 2.0, sp):channel(e0 / 1.5, sp)
    return [mean(sr[rr]) / mean(sr[ss]) - 1.0 for sr in srs]
end

"""
    multirank(specs::AbstractArray{<:Spectrum})::Float64

A single number that compares the low and high energy portions of the spectra for similarity.  A `multirank(...)` score
of zero means all the spectra are very similar and a large number means very different. A high score suggests that
one or more of the spectra may suffer from additional low energy absorption due to surface roughness, an obstruction,
sample tilt or other.
"""
function multirank(specs::AbstractArray{<:Spectrum})::Float64
    -(-)(extrema(multiscore(specs))...)
end

# Return the common prefix
function _commonname(nms)
    alleq(cs) = all(c->c==cs[1], cs[2:end])
    l = minimum(length, nms)
    for i in Base.OneTo(l)
        if !alleq(map(n->n[nextind(n, 0, i)], nms))
            return i>1 ? nms[1][1:nextind(nms[1],0, i-1)] : "Multisum[$(spectrumCounter())["
        end
    end
    return nms[1][1:l]*"["
end

"""
    multisum(specs::Spectrum{T})::Spectrum{T} where {T <: Real}

Sum together spectra collected from multiple detectors simultaneously from
the same electrons interacting with the same material for the real-time.  
The detectors should be calibrated close to identically to maintain
the detector resolution and peak positions.
"""
function multisum(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where {T <: Real}
    checkmultispec(specs)
    sc = [ sum(ss) for ss in zip(map(s->s.counts, specs)...) ]
    cp = commonproperties(specs)
    cp[:LiveTime] = sum(s->s[:LiveTime], specs)
    cp[:Name] = _commonname(map(s->s[:Name], specs))*"sum]"
    return Spectrum(specs[1].energy, sc, cp)
end

"""
    multimean(specs::Spectrum{T})::Spectrum{T} where {T <: Real}

Average on a channel-by-channel basis spectra collected from multiple detectors
simultaneously from the same electrons interacting with the same material for 
the real-time.  The detectors should be calibrated close to identically to maintain
the detector resolution and peak positions.
"""
function multimean(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where {T <: Real}
    checkmultispec(specs)
    sc = [ mean(ss) for ss in zip(map(s->s.counts, specs)...) ]
    cp = commonproperties(specs)
    cp[:LiveTime] = sum(s->s[:LiveTime], specs)
    cp[:Name] = _commonname(map(s->s[:Name], specs))*"mean]"
    return Spectrum(specs[1].energy, sc, cp)
end

"""
    multicompare(specs::AbstractArray{Spectrum{T}}) where {T <: Real}

Compares the intensity for the spectra in `specs` against the mean intensity
on a channel-by-channel basis.  Compute the ratio for each channel in each 
spectrum of the spectrum intensity over the mean intensity for that channel.
You expect the ratio to be unity when the spectra are identical and deviate
from unity when the spectra are different. 
"""
function multicompare(specs::AbstractArray{Spectrum{T}}) where {T <: Real}
    checkmultispec(specs)
    hs=sqs = map(ss -> sqrt.(max.(ss.counts, 1)), specs)
    sc = [ mean(ss) for ss in zip(sqs...) ]
    return map(sq -> (sum(sc)/sum(sq))*(sq./sc), sqs)
end