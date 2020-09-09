
struct FilterFitPacket
    detector::Detector
    filter::TopHatFilter
    references::Vector{FilteredReference}
end # struct

function Base.show(io::IO, ffp::FilterFitPacket)
    lines = join(map(r->"\t$(r.identifier),", ffp.references),"\n")
    print(io, "References[\n\t$(ffp.detector), \n$lines\n]")
end

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

references(refs::AbstractVector{<:Tuple{Spectrum, Element, Material}}, fwhm::Float64)::FilterFitPacket =
    references(refs, matching(refs[1][1], fwhm))

function references(refs::AbstractVector{<:Tuple{Spectrum, Element, Material}}, det::EDSDetector)::FilterFitPacket
    chcount = length(refs[1][1])
    @assert all(length(r[1])==chcount for r in refs[2:end]) "The number of channels must match in all the spectra."
    @assert all(matches(r[1], det) for r in refs[1:end]) "All the spectra must match the detector $det."
    ff = buildfilter(det)
    refs = mapreduce(ref->filterreference(ff, ref...), append!, refs)
    return FilterFitPacket(det, ff, selectBestReferences(refs))
end

fit(spec::Spectrum, ffp::FilterFitPacket) = fit(FilteredUnknownW, spec, ffp.filter, ffp.refs)
"""
fit(hs::HyperSpectrum, ffp::FilterFitPacket; mode::Symbol=:Fast, zero = x -> max(0.0, x))::Array{KRatios}

  * mode = :Fast - Uses precomputed, filtered "vector" fit method
  * mode = :Intermediate - Uses an optimized full fit without refits for k(s) < 0.0
  * mode = :Full - Uses the full single spectrum fit algorithm including refitting when one or more k < 0.0

Performs a filtered fit on a hyperspectrum returning an `Array{KRatios}`.

Selecting a mode:
  :Fast is good for generating k-ratio maps or exploratory analysis of a k-ratio map. :Full is best when a
  quantitative map of a high count hyperspectrum is desired.  Fit results for major elements are similar for
  all three but differ for minor and trace elements.  Particularly when a k-ratio is slightly negative. This
  negative k-ratio can effect the other k-ratios.  :Fast also works less well when many elements (>>10) (particularly
  interfering elements) are included in the fit. Unfortunately, :Intermediate and :Full slow down when many elements
  are fit - O(n(elements)²).

The following timing on a 512 x 512 x 2048 hyperspectrum fitting 15 elements with 25 distinct ROIs on a fast laptop
with 64 GiB memory give a relative feel for the speed of each algorithm.  Yes, :Fast is approximately 20x faster than
:Intermediate and used almost 100x less memory.  (Single thread timings)

|---------------|----------------|--------------|--------------|----------|
| mode          | Run time (s)   | Allocations  | Memory (GiB) | GC time  |
|---------------|----------------|--------------|--------------|----------|
| :Fast         | 44.8           | 3.95 M       | 5.72         | 4.2%     |
| :Intermediate | 1305.1         | 24.83 M      | 523.7        | 2.7%     |
| :Full         | 2056           | 2.7 G        | 786.1        | 4.2%     |
|---------------|----------------|--------------|--------------|----------|
"""
function fit(hs::HyperSpectrum, ffp::FilterFitPacket; mode::Symbol=:Fast, zero = x -> max(0.0, x))::Array{KRatios}
    _tophatfilterhs = (hs, data, thf, scale) -> begin
        @assert length(data)<=length(thf) "The reference spectra have fewer channels than the hyperspectrum data."
        filtered = Float64[ filtereddatum(thf,data,i) for i in eachindex(data) ]
        dp = Float64[ max(x, 1.0) for x in data ] # To ensure covariance isn't zero or infinite precision
        covar = Float64[ filteredcovar(thf, dp, i, i) for i in eachindex(data) ]
        return FilteredUnknownW(UnknownLabel(hs), scale, eachindex(data), data, filtered, covar)
    end
    fitcontiguousx = (unk, ffs, chs) ->
        _buildscale(unk, ffs) * pinv(_buildmodel(ffs, chs), rtol=1.0e-6)*extract(unk, chs)
    _filterfitx = (unk, ffs, fitrois) ->
        Iterators.flatten(map(fr -> fitcontiguousx(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs), fr), fitrois))
    @assert matches(hs[CartesianIndices(hs)[1]], ffp.detector) "The detector for the hyper-spectrum must match the detector for the filtered references."
    if mode==:Fast
        vq = VectorQuant(ffp.references, ffp.filter)
        return fit(vq, hs, zero)
    elseif mode==:Intermediate
        krs = zeros(Float32, length(ffp.references), size(hs)...)
        idose = 1.0/dose(hs)
        len = 1:depth(hs)
        data = zeros(Float64, length(ffp.filter))
        fitrois = ascontiguous(map(fd -> fd.ffroi, ffp.references))
        for ci in CartesianIndices(hs)
            data[len] = hs.counts[len, ci]
            unk = _tophatfilterhs(hs, data, ffp.filter, idose)
            krs[:, ci] .= zero.(_filterfitx(unk, ffp.references, fitrois))
        end
        res = KRatios[]
        for i in filter(ii->ffp.references[ii].identifier isa CharXRayLabel, eachindex(ffp.references))
            k, lbl = krs[i], ffp.references[i].identifier
            rprops = properties(spectrum(lbl))
            push!(res, KRatios(xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i,:,:]))
        end
        return res
    elseif mode==:Full
        krs = zeros(Float32, length(ffp.references), size(hs)...)
        idose = 1.0/dose(hs)
        len = 1:depth(hs)
        data = zeros(Float64, length(ffp.filter))
        for ci in CartesianIndices(hs)
            data[len] = hs.counts[len, ci]
            unk = _tophatfilterhs(hs, data, ffp.filter, idose)
            uvs = _filterfit(unk, ffp.references)
            krs[:, ci] = map(id->convert(Float32, value(id, uvs)), ( ref.identifier for ref in ffp.references ))
        end
        res = KRatios[]
        for i in filter(ii->ffp.references[ii].identifier isa CharXRayLabel, eachindex(ffp.references))
            k, lbl = krs[i], ffp.references[i].identifier
            rprops = properties(spectrum(lbl))
            push!(res, KRatios(xrays(lbl), properties(hs), rprops, rprops[:Composition], krs[i,:,:]))
        end
        return res
    else
        @assert false "The mode argument must be :Fast, :Intermediate or :Full."
    end
end
