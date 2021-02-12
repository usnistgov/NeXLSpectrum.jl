# Routines to handle spectra collected simulataneously on multiple detectors

using Formatting

"""
    loadmultispec(path::AbstractString, basefn::AbstractString, indexes=0:3, fnmapper::String = "{1}[{2}].msa")

Load multiple spectra using the `basefn` and `fnmapper` to determine which spectra to load.  The spectra should be
related in the sense that they were all collected simulataneously so they have the same `:RealTime`, `:BeamEnergy` 
and `:LiveTime`.

"""
function loadmultispec(
    path::AbstractString,
    basefn::AbstractString,
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
    @assert all(
        get(sp, :BeamEnergy, 0.0) == get(specs[1], :BeamEnergy, 0.0) for sp in specs[2:end]
    ) "All the beam energies must match."
    @assert all(
        get(sp, :RealTime, 0.0) == get(specs[1], :RealTime, 0.0) for sp in specs[2:end]
    ) "All the real times must match."
    @assert all(
        get(sp, :ProbeCurrent, 0.0) == get(specs[1], :ProbeCurrent, 0.0) for
        sp in specs[2:end]
    ) "All the probe currents must match."
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