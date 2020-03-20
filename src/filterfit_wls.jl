using SparseArrays
using Polynomials
using LinearAlgebra

"""
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.  This is an approximation that produces over
optimistic resulting covariance matrix.
"""
struct FilteredUnknownW <: FilteredUnknown
    identifier::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data (always 1:...)
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::Vector{Float64} # Channel covariance
end

# _buildXXX - Helpers for fitcontiguousX functions...
function _buildmodel(ffs::Array{FilteredReference}, chs::UnitRange{Int})::Matrix{Float64}
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    return x
end

"""
    covariance(fd::FilteredUnknownW, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance diagnonal elements over the specified range of channels.
<code>roi</code> must be fully contained within the data in <code>fd</code>.
"""
NeXLUncertainties.covariance(fd::FilteredUnknownW, roi::UnitRange{Int})::AbstractVector{Float64} =
    fd.covariance[roi]

"""
Weighted least squares for FilteredUnknownW
"""
function fitcontiguousww(unk::FilteredUnknownW, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    wgts = [ff.covscale for ff in ffs]
    return scale * wlspinv(extract(unk, chs), x, covariance(unk, chs), wgts, lbls)
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
    filter(spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknown

For filtering the unknown spectrum. Defaults to the weighted fitting model.
"""
Base.filter(spec::Spectrum, filt::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknown =
    filter(FilteredUnknownW, spec, filt, scale)

"""
    filter(::Type{FilteredUnknownW}, spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknownW

For filtering the unknown spectrum. Process the full Spectrum with the specified filter for use with the weighted
least squares model.
"""
function Base.filter(::Type{FilteredUnknownW}, spec::Spectrum, filter::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownW
    range(i) = filter.offsets[i]:filter.offsets[i]+length(filter.filters[i])-1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [dot(filter.filters[i], view(data, range(i))) for i in eachindex(data)]
    covariance(filter, data) = # The diagnonal elements of the full covariance matrix
        [dot(filter.filters[r] .* view(data, range(r)), filter.filters[r]) for r in eachindex(data)]
    data = counts(spec, Float64, true)
    # Compute the filtered data
    filtered = apply(filter, data)
    roi = eachindex(filtered)
    # max(d,1.0) is necessary to ensure the variances are positive
    covar = covariance(filter, map(d -> max(d, 1.0), data))
    return FilteredUnknownW(UnknownLabel(spec), scale, roi, roi, data, filtered, covar)
end

"""
    filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg=fitcontiguousww)::UncertainValues

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
function filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg = fitcontiguousww, forcezeros = true)::FilterFitResult
    trimmed, refit, removed, retained = copy(ffs), true, Vector{UncertainValues}(), nothing # start with all the FilteredReference
    while refit
        refit = false
        fitrois = ascontiguous(map(fd -> fd.ffroi, trimmed))
        retained = map(fr -> alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr), fitrois)
        kr = cat(retained)
        if forcezeros
            for lbl in keys(kr)
                if value(lbl, kr) < 0.0
                    splice!(trimmed, findfirst(ff -> ff.identifier == lbl, trimmed))
                    push!(removed, uvs([lbl], [0.0], reshape([σ(lbl, kr)], (1, 1))))
                    refit = true
                end
            end
        end
    end # while
    kr = cat(append!(retained, removed))
    resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, resid, pb)
end

fit(ty::Type{FilteredUnknownW}, unk::Spectrum, filt::TopHatFilter, refs::Vector{FilteredReference}, forcezeros = true) =
    filterfit(filter(ty, unk, filt, 1.0 / dose(unk)), refs, fitcontiguousww, forcezeros)

fit(unk::Spectrum, filt::TopHatFilter, refs::Vector{FilteredReference}, forcezeros = true) =
    fit(FilteredUnknownW, unk, filt, refs, forcezeros)
