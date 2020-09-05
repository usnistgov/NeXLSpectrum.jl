


struct FilterFitPacket
    detector::Detector
    filter::TopHatFilter
    references::Vector{FilteredReference}
end # struct

function reference(elm::Element, spec::Spectrum, mat::Material; pc = nothing, lt = nothing, e0 = nothing)
    if !isnothing(lt)
        @assert lt>0.0 "The live time must be larger than zero."
        spec[:LiveTime] = lt
    end
    if !isnothing(pc)
        @assert pc>0.0 "The probe current must be larger than zero."
        spec[:ProbeCurrent] = pc
    end
    if !isnothing(e0)
        @assert e0 >= 1.0e3 "The beam energy ($(e0/1000.0)) keV is too low."
        spec[:BeamEnergy] = e0
    end
    @assert haskey(spec, :LiveTime) "The :LiveTime property must be defined.  (Use the `lt` keyword argument)"
    @assert haskey(spec, :ProbeCurrent) "The :ProbeCurrent property must be defined.  (Use the `pc` keyword argument)"
    @assert haskey(spec, :BeamEnergy) "The :BeamEnergy property must be defined.  (Use the `e0` keyword argument)"
    return ( spec, elm, mat)
end

reference(elm::Element, spec::Spectrum) = reference(elm, spec, spec[:Composition]; kwargs)

reference(elm::Element, filename::AbstractString, mat::Material; kwargs...) =
    reference( elm, loadspectrum(filename), mat; kwargs... )

reference(elm::Element, filename::AbstractString; kwargs...) =
    reference( elm, loadspectrum(filename); kwargs... )

function references(refs::AbstractVector{<:Tuple{Spectrum, Element, Material}}, fwhm::Float64)::FilterFitPacket
    chcount = length(refs[1][1])
    @assert all(length(r[1])==chcount for r in refs[2:end]) "The number of channels must match in all the spectra."
    det = matching(refs[1][1], fwhm)
    @assert all(matches(r[1], det) for r in refs[2:end]) "All the spectra must match the detector $det."
    ff = buildfilter(det)
    return FilterFitPacket(det, ff, mapreduce(ref->filterreference(ff, ref...), append!, refs))
    filterreference()
end

fit(spec::Spectrum, ffp::FilterFitPacket) = fit(FilteredUnknownW, spec, ffp.filter, ffp.refs)

function fit(hs::HyperSpectrum, ffp::FilterFitPacket; fast=true, zero = x -> max(0.0, x))::Array{KRatios}
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector)
    if fast
        vq = VectorQuant(ffp.references, ffp.filter)
        return fit(vq, hs, zero)
    else
        @assert false "Not yet implemented."
    end
end
