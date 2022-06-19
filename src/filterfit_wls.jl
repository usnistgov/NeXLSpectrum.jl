
"""
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.  This is an approximation that produces over
optimistic resulting covariance matrix.
"""
struct FilteredUnknownW{T} <: FilteredUnknown where { T<:AbstractFloat}
    label::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data (always 1:...)
    data::Vector{T} # Spectrum data over ffroi
    filtered::Vector{T} # Filtered spectrum data over ffroi
    covariance::Vector{T} # Channel covariance
end

# _buildXXX - Helpers for fitcontiguousX functions...
function _buildmodel(
    ffs::AbstractVector{FilteredReference{T}},
    roi::UnitRange{Int},
) where { T<: AbstractFloat }
    x = zeros(T, length(roi), length(ffs))
    data = zeros(T, length(roi))
    for (i,fd) in enumerate(ffs)
        #@assert fd.ffroi.start >= roi.start "$(fd.ffroi.start) < $(roi.start) in $(fd)"
        #@assert fd.ffroi.stop <= roi.stop "$(fd.ffroi.stop) > $(roi.stop) in $(fd)"
        fill!(data, zero(T))
        data[fd.ffroi.start-roi.start+1:fd.ffroi.stop-roi.start+1] .= fd.filtered
        x[:, i] .= data
        # @assert x[:, i] == extract(fd, roi)
    end
    return x
end

"""
Weighted least squares for FilteredUnknownW
"""
function fitcontiguousww(
    unk::FilteredUnknownW{T},
    ffs::AbstractVector{FilteredReference{T}},
    chs::UnitRange{Int},
)::UncertainValues where { T <: AbstractFloat }
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    # dcs is a factor that accounts for the heteroskedasciscity introduced by the filter
    dcs = Diagonal([ T(ff.covscale) for ff in ffs ])
    w = Diagonal([sqrt(one(T) / T(cv)) for cv in view(unk.covariance, chs)])
    genInv = pinv(w * x, rtol = 1.0e-6)
    return scale * uvs(lbls, genInv * w * extract(unk, chs), dcs * (genInv * transpose(genInv)) * dcs)
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
tophatfilter(spec::Spectrum, filt::TopHatFilter{T}, scale::Float64 = 1.0) where { T <: AbstractFloat } = #
    tophatfilter(FilteredUnknownW{T}, spec, filt, scale)

"""
    tophatfilter(::Type{FilteredUnknownW}, spec::Spectrum, thf::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknownW

For filtering the unknown spectrum. Process the full Spectrum with the specified filter for use with the weighted
least squares model.
"""
function tophatfilter(
    ::Type{FilteredUnknownW{T}},
    spec::Spectrum,
    thf::TopHatFilter{T},
    scale::Float64 = 1.0,
) where { T<: AbstractFloat }
    data = counts(spec, 1:length(thf), T, true)
    filtered = T[filtereddatum(thf, data, i) for i in eachindex(data)]
    dp = T[max(x, 1.0) for x in data] # To ensure covariance isn't zero or infinite precision
    covar = T[filteredcovar(thf, dp, i, i) for i in eachindex(data)]
    return FilteredUnknownW{T}(
        UnknownLabel(spec),
        scale,
        1:length(data),
        data,
        filtered,
        covar,
    )
end

# Actually perform the filter fit and return k-ratios in an UncertainValues object
function _filterfit(
    unk::FilteredUnknownW{T},
    ffs::AbstractVector{FilteredReference{T}},
    forcezeros,
) where { T <: AbstractFloat }
    trimmed, refit, removed = copy(ffs), true, UncertainValues[] # start with all the FilteredReference
    while refit
        refit = false
        fitrois = ascontiguous(map(fd -> fd.ffroi, trimmed))
        retained = map(fitrois) do fr
            # `fitcontiguousww(..) performs the fit
            fitcontiguousww(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr)
        end
        kr = cat(retained)
        if forcezeros
            for lbl in keys(kr)
                if NeXLUncertainties.value(kr, lbl) < 0.0
                    splice!(trimmed, findfirst(ff -> ff.label == lbl, trimmed))
                    push!(removed, uvs([lbl], [0.0], reshape([σ(kr, lbl)], (1, 1))))
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
function filterfit(
    unk::FilteredUnknownW{T},
    ffs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
)::FilterFitResult where { T <: AbstractFloat }
    krs = _filterfit(unk, ffs, forcezeros)
    resid, pb = _computeResidual(unk, ffs, krs), _computecounts(unk, ffs, krs)
    return FilterFitResult(unk.label, krs, unk.roi, unk.data, resid, pb)
end

function fit_spectrum(
    ty::Type{FilteredUnknownW{T}},
    unk::Spectrum,
    filt::TopHatFilter,
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat }
    bestRefs = selectBestReferences(refs)
    return filterfit(tophatfilter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, forcezeros)
end

function fit_spectrum(
    ty::Type{FilteredUnknownW{T}},
    unks::AbstractVector{<:Spectrum},
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat }
    bestRefs = selectBestReferences(refs)
    return map(unks) do unk
        fu = tophatfilter(ty, unk, filt, 1.0 / dose(unk))
        filterfit(fu, bestRefs, forcezeros)
    end
end

fit_spectrum(
    unk::Spectrum,
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat } = #
    fit_spectrum(FilteredUnknownW{T}, unk, filt, refs, forcezeros)

fit_spectrum(
    unks::AbstractVector{Spectrum},
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat } = # 
    fit_spectrum(FilteredUnknownW{T}, unks, filt, refs, forcezeros)
