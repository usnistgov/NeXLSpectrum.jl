"""
    FilteredUnknownG

Represents the unknown in a filter fit using the full generalized fitting model.  This model is expensive to
calculate but uses the full generalized linear fitting model which produces the correct fit covariances.
"""
struct FilteredUnknownG <: FilteredUnknown
    identifier::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data (always 1:...)
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::AbstractMatrix{Float64} # Channel covariance
end

"""
    filter(spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredDatum

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
function Base.filter(::Type{FilteredUnknownG}, spec::Spectrum, filter::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownG
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
    function covariance(roi, data, filter) # Calculated the full covariance matrix
        off(r, o) = r.start-o+1:r.stop-o+1
        rs, cs, vs = Int[], Int[], Float64[]
        res = zeros( ( length(roi), length(roi) ) )
        bdata = map(d -> max(d, 1.0), data)
        for r in roi
            rr, rf, ro = range(r), filter.filters[r], filter.offsets[r]
            for c in roi
                ii = intersect(rr, range(c))
                if length(ii) > 0
                    cf, co = filter.filters[c], filter.offsets[c]
                    # The next line is the core time sink...
                    v = dot( view(rf, off(ii, ro)) .* view(bdata, ii) , view(cf, off(ii, co)))
                    if v ≠ 0.0
                        push!(rs, r), push!(cs, c), push!(vs, v)
                    end
                end
            end
        end
        push!(rs, roi.stop), push!(cs, roi.stop), push!(vs, 0.0) # force the size
        return sparse(rs, cs, vs)
    end
    data = counts(spec, Float64, true)
    # Compute the filtered data
    filtered = apply(filter, data)
    roi = 1:findlast(f -> f ≠ 0.0, filtered)
    # max(d,1.0) is necessary to ensure the variances are positive
    # diag = sparse(roi,roi,map(d -> max(d, 1.0), data[roi]))
    covar = covariance(roi, data, filter)
    @assert covar isa AbstractSparseMatrix "The covariance matrix should be a sparse matrix in filter(...)::FilterUnknown."
    checkcovariance!(covar)
    return FilteredUnknownG(UnknownLabel(spec), scale, roi, roi, data[roi], filtered[roi], covar)
end

"""
    covariance(fd::FilteredUnknownG, roi::UnitRange{Int})::AbstractMatrix{Float64}

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.  <code>roi</code> must
be fully contained within the filtered edata in <code>fd</code>.
"""
NeXLUncertainties.covariance(fd::FilteredUnknownG, roi::UnitRange{Int})::AbstractMatrix{Float64} =
     fd.covariance[roi,roi]

"""
    filterfit(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, alg = fitcontiguousp, forcezeros::Bool=true)::FilterFitResult
 trimmed, refit, removed, retained = copy(ffs), true, UncertainValues[], nothing # start with all the FilteredReference
 while refit
     fitrois = ascontiguous(map(fd->fd.ffroi, trimmed))
     # @info "Fitting $(length(trimmed)) references in $(length(fitrois)) blocks - $fitrois"
     retained = map(fr->alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr), fitrois)
     kr = cat(retained)
     refit = false
     for lbl in keys(kr)
         if value(lbl, kr) <= 0.0
             splice!(trimmed, findfirst(ff->ff.identifier==lbl, trimmed))
             push!(removed, uvs([lbl],[forcezeros ? 0.0 : value(lbl, kr)], reshape([σ(lbl, kr)], (1,1))))
             refit=true
         end
     end
 end # while
 kr = cat(append!(retained, removed))
 resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
 return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, resid, pb)
end

"""
Generalized least squares (my implementation)
"""
function fitcontiguousg(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glssvd(extract(unk, chs), x, covariance(unk, chs), lbls)
end

"""
Generalized least squares (pseudo-inverse)
"""
function fitcontiguousp(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, Matrix(covariance(unk, chs)), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, Diagonal(diag(covariance(unk, chs))), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw2(unk::FilteredUnknownG, ffs::AbstractVector{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * wlssvd(extract(unk, chs), x, diag(covariance(unk, chs)), lbls)
end

function fit(ty::Type{FilteredUnknownG}, unk::Spectrum, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true)
    bestRefs = selectBestReferences(refs)
    return filterfit(filter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, fitcontiguousww, forcezeros)
end

function fit(ty::Type{FilteredUnknownG}, unks::AbstractVector{Spectrum}, filt::TopHatFilter, refs::AbstractVector{FilteredReference}, forcezeros = true)
    bestRefs = selectBestReferences(refs)
    return map(unk->filterfit(filter(ty, unk, filt, 1.0 / dose(unk)), bestRefs, fitcontiguousww, forcezeros), unks)
end
