"""
    FilteredUnknownG

Represents the unknown in a filter fit using the full generalized fitting model.  This model is expensive to
calculate but uses the full generalized linear fitting model which produces the correct fit covariances.
"""
struct FilteredUnknownG <: FilteredUnknown
    label::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    data::Vector{Float64} # Spectrum data
    filtered::Vector{Float64} # Filtered spectrum data
    covariance::AbstractMatrix{Float64} # Channel covariance
end

"""
    tophatfilter(spec::Spectrum, thf::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredDatum

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
function tophatfilter(::Type{FilteredUnknownG}, spec::Spectrum, thf::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownG
    data = counts(spec, Float64, true)
    filtered = [ filtereddatum(thf, data, i) for i in eachindex(data) ]
    dp = map(x->max(x, 1.0), data) # To ensure covariance isn't zero or infinite precision
    covar = [ filteredcovar(thf, dp, r, c) for r in eachindex(data), c in eachindex(data) ]
    return FilteredUnknownG(UnknownLabel(spec), scale, data, filtered, covar)
end

"""
    covariance(fd::FilteredUnknownG, roi::UnitRange{Int})::AbstractMatrix{Float64}

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.  `roi` must
be fully contained within the filtered edata in `fd`.
"""
NeXLUncertainties.covariance(fd::FilteredUnknownG, roi::UnitRange{Int})::AbstractMatrix{Float64} =
    fd.covariance[roi, roi]

"""
    filterfit(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, alg::Function = fitcontiguousp, forcezeros::Bool=true)::FilterFitResult
    trimmed, removed = copy(ffs), UncertainValues[] # start with all the FilteredReference
    while true
        # Build a list of contiguous rois each of which is fit as a batch
        fitrois = ascontiguous(map(fd->fd.ffroi, trimmed))
        # @info "Fitting $(length(trimmed)) references in $(length(fitrois)) blocks - $fitrois"
        retained = map(fitrois) do fr
            alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr)
        end
        kr = cat(retained)
        refit = false
        for lbl in keys(kr)
            if NeXLUncertainties.value(lbl, kr) <= 0.0
                splice!(trimmed, findfirst(ff->ff.label==lbl, trimmed))
                push!(removed, uvs([lbl],[forcezeros ? 0.0 : NeXLUncertainties.value(lbl, kr)], reshape([σ(lbl, kr)], (1,1))))
                refit=true
            end
        end
        if !refit
            kr = cat(append!(retained, removed))
            resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
            return FilterFitResult(unk.label, kr, eachindex(unk.data), unk.data, resid, pb)        
        end
    end
end

"""
Generalized least squares (my implementation)
"""
function fitcontiguousg(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glssvd(extract(unk, chs), x, unk.covariance[chs,chs], lbls)
end

"""
Generalized least squares (pseudo-inverse)
"""
function fitcontiguousp(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, unk.covariance[chs,chs], lbls)
end

"""
Generalized least squares (inverse)
"""
function fitcontiguousi(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glsinv(extract(unk, chs), x, unk.covariance[chs,chs], lbls)
end

"""
Generalized least squares (Cholesky)
"""
function fitcontiguousc(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glschol(extract(unk, chs), x, unk.covariance[chs,chs], lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance as a matrix fed into glspinv.
"""
function fitcontiguousw(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, Diagonal(diag(unk.covariance[chs,chs])), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance as a vector fed into wlssvd.
"""
function fitcontiguousw2(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * wlssvd(extract(unk, chs), x, diag(unk.covariance[chs,chs]), lbls)
end

function fit_spectrum(ty::Type{FilteredUnknownG}, unk::Spectrum, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros::Bool = true)
    bestRefs = selectBestReferences(refs)
    return filterfit(tophatfilter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, fitcontiguousc, forcezeros)
end

function fit_spectrum(ty::Type{FilteredUnknownG}, unks::AbstractVector{<:Spectrum}, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros::Bool = true)
    bestRefs = selectBestReferences(refs)
    return map(unk->filterfit(tophatfilter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, fitcontiguousc, forcezeros), unks)
end


