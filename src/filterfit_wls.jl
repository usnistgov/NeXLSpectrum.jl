
"""
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.  This is an approximation that produces over
optimistic resulting covariance matrix.
"""
struct FilteredUnknownW <: FilteredUnknown
    label::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data (always 1:...)
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::Vector{Float64} # Channel covariance
end

# _buildXXX - Helpers for fitcontiguousX functions...
function _buildmodel(ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::Matrix{Float64}
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    return x
end

"""
    covariance(fd::FilteredUnknownW, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance diagnonal elements over the specified range of channels.
`roi` must be fully contained within the data in `fd`.
"""
NeXLUncertainties.covariance(fd::FilteredUnknownW, roi::UnitRange{Int})::AbstractVector{Float64} =
    view(fd.covariance,roi)

"""
Weighted least squares for FilteredUnknownW
"""
function fitcontiguousww(unk::FilteredUnknownW, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    covscales = [ff.covscale for ff in ffs]
    return scale * wlspinv2(extract(unk, chs), x, covariance(unk, chs), covscales, lbls)
end

function ascontiguous(rois::AbstractArray{UnitRange{Int}})
    # Join the UnitRanges into contiguous UnitRanges
    join(roi1, roi2) = min(roi1.start, roi2.start):max(roi1.stop, roi2.stop)
    srois = sort(rois)
    res = [srois[1]]
    for roi in srois[2:end]
        if length(intersect(res[end], roi)) > 0
            res[end] = join(roi, res[end]) # Join UnitRanges
        else
            push!(res, roi) # Add a new UnitRange
        end
    end
    return res
end

"""
    tophatfilter(spec::Spectrum, thf::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknown

For filtering the unknown spectrum. Defaults to the weighted fitting model.
"""
tophatfilter(spec::Spectrum, filt::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknown =
    tophatfilter(FilteredUnknownW, spec, filt, scale)

"""
    tophatfilter(::Type{FilteredUnknownW}, spec::Spectrum, thf::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknownW

For filtering the unknown spectrum. Process the full Spectrum with the specified filter for use with the weighted
least squares model.
"""
function tophatfilter(::Type{FilteredUnknownW}, spec::Spectrum, thf::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownW
    data = counts(spec, 1:length(thf), Float64, true)
    filtered = Float64[ filtereddatum(thf,data,i) for i in eachindex(data) ]
    dp = Float64[ max(x, 1.0) for x in data ] # To ensure covariance isn't zero or infinite precision
    covar = Float64[ filteredcovar(thf, dp, i, i) for i in eachindex(data) ]
    return FilteredUnknownW(UnknownLabel(spec), scale, eachindex(data), data, filtered, covar)
end

# Actually perform the filter fit and return k-ratios in an UncertainValues object
function _filterfit(unk::FilteredUnknownW, ffs::AbstractVector{FilteredReference}, forcezeros)
    trimmed, refit, removed = copy(ffs), true, UncertainValues[] # start with all the FilteredReference
    while refit
        refit = false
        fitrois = ascontiguous(map(fd -> fd.ffroi, trimmed))
        # `fitcontiguousww(..) performs the fit
        retained = map(fr -> fitcontiguousww(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr), fitrois)
        kr = cat(retained)
        if forcezeros
            for lbl in keys(kr)
                if NeXLUncertainties.value(lbl, kr) < 0.0
                    splice!(trimmed, findfirst(ff -> ff.label == lbl, trimmed))
                    push!(removed, uvs([lbl], [0.0], reshape([σ(lbl, kr)], (1, 1))))
                    refit = true
                end
            end
        end
        if !refit
            return cat(append!(retained, removed))
        end
    end # while
    @assert false
    return removed # To maintain type
end

"""
    filterfit(unk::FilteredUnknownW, ffs::AbstractVector{FilteredReference}, forcezeros = true)::FilterFitResult

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).

This function is designed to reperform the fit if one or more k-ratio is less-than-or-equal-to zero.  The
FilteredReference corresponding to the negative value is removed from the fit and the fit is reperformed. How the
non-positive value is handled is determine by forcezeros. If forcezeros=true, then the returned k-ratio for the
non-positive value will be set to zero (but the uncertainty remains the fitted one).  However, if forcezeros=false,
then the final non-positive k-ratio is returned along with the associated uncertainty.  forcezeros=false is better
when a number of fit k-ratio sets are combined to produce an averaged k-ratio with reduced uncertainty. forcezeros=true
would bias the result positive.
"""
function filterfit(unk::FilteredUnknownW, ffs::AbstractVector{FilteredReference}, forcezeros = true)::FilterFitResult
    kr = _filterfit(unk, ffs, forcezeros)
    resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
    return FilterFitResult(unk.label, kr, unk.roi, unk.data, resid, pb)
end

function filterfitk(unk::FilteredUnknownW, ffs::AbstractVector{FilteredReference}, forcezeros = true)
    uvs = _filterfit(unk, ffs, forcezeros)
    unkprops = properties(spectrum(unk.label))
    return [ KRatio(xrays(lbl), unkprops, properties(spectrum(lbl)), spectrum(lbl)[:Composition], uvs[lbl])
        for lbl in labels(uvs) ]
end

function fit_spectrum(ty::Type{FilteredUnknownW}, unk::Spectrum, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true)
    bestRefs = selectBestReferences(refs)
    return filterfit(tophatfilter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, forcezeros)
end

function fit_spectrum(ty::Type{FilteredUnknownW}, unks::AbstractVector{<:Spectrum}, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true)
    bestRefs = selectBestReferences(refs)
    return map(unk->filterfit(tophatfilter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, forcezeros), unks)
end

fit_spectrum(unk::Spectrum, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true) =
    fit_spectrum(FilteredUnknownW, unk, filt, refs, forcezeros)

fit_spectrum(unks::AbstractVector{Spectrum}, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true) =
    fit_spectrum(FilteredUnknownW, unks, filt, refs, forcezeros)
