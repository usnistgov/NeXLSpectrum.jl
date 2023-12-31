
"""
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.  This is an approximation that produces over
optimistic resulting covariance matrix.
"""
struct FilteredUnknownW{T} <: FilteredUnknown where { T<:AbstractFloat}
    label::UnknownLabel # A way of identifying this filtered datum
    scale::T # A dose or other scale correction factor
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
    for (i,fd) in enumerate(ffs)
        #@assert fd.ffroi.start >= roi.start "$(fd.ffroi.start) < $(roi.start) in $(fd)"
        #@assert fd.ffroi.stop <= roi.stop "$(fd.ffroi.stop) > $(roi.stop) in $(fd)"
        x[fd.ffroi.start-roi.start+1:fd.ffroi.stop-roi.start+1, i] .= fd.filtered
        # @assert x[:, i] == extract(fd, roi)
    end
    return x
end

"""
    fitcontiguousww(
        unk::FilteredUnknownW{T},
        ffs::AbstractVector{FilteredReference{T}},
        chs::UnitRange{Int},
    )::UncertainValues

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
    genInv = pinv(w * x, rtol = 1.0e-6) # 60% of allocation here
    ext = extract(unk, chs)
    # @assert all(eltype.((genInv, w, dcs, ext)).==T)
    # @assert eltype(scale)==Float64
    return scale * uvs(lbls, genInv * w * ext, dcs * (genInv * transpose(genInv)) * dcs)
end

"""
    fitcontiguouso(
        unk::FilteredUnknown,
        ffs::AbstractVector{FilteredReference{T}},
        chs::UnitRange{Int},
    )::UncertainValues
    
Ordinary least squares for either FilteredUnknown[G|W]
"""
function fitcontiguouso(
    unk::FilteredUnknown,
    ffs::AbstractVector{FilteredReference{T}},
    chs::UnitRange{Int},
)::UncertainValues where { T <: AbstractFloat }
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    genInv = pinv(x, rtol = 1.0e-6)
    return scale * uvs(lbls, genInv * extract(unk, chs), genInv * transpose(genInv))
end

function ascontiguous(rois::AbstractArray{UnitRange{Int}})::Vector{UnitRange{Int}}
    # Join the UnitRanges into contiguous UnitRanges
    join(roi1, roi2) = min(roi1.start, roi2.start):max(roi1.stop, roi2.stop)
    srois = sort(rois)
    res = UnitRange{Int}[srois[1]]
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
    tophatfilter(spec::Spectrum, filt::TopHatFilter{T}, scale::T = one(T))::FilteredUnknown where { T <: AbstractFloat }

For filtering the unknown spectrum. Defaults to the weighted fitting model.
"""
tophatfilter(spec::Spectrum, filt::TopHatFilter{T}, scale::T = one(T)) where { T <: AbstractFloat } = #
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
    scale::T = one(T),
) where { T<: AbstractFloat }
    data = counts(spec, 1:length(thf), T, true)
    filtered = T[filtereddatum(thf, data, i) for i in eachindex(data)]
    dp = T[max(x, 1.0) for x in data] # To ensure covariance isn't zero or infinite precision
    covar = T[filteredcovar(thf, dp, i, i) for i in eachindex(data)]
    return FilteredUnknownW{T}(
        UnknownLabel(spec),
        T(scale),
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
    forcezeros::Bool,
)::UncertainValues where { T <: AbstractFloat }
    trimmed, removed = copy(ffs), UncertainValues[] # start with all the FilteredReference
    while true
        retained = map(ascontiguous(map(fd -> fd.ffroi, trimmed))) do fr
            # `fitcontiguousww(..) performs the fit
            fitcontiguousww(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr)
        end
        if forcezeros
            vals  = Dict{Label, UncertainValue}()
            for kr in retained, lbl in filter(l->value(kr,l) < 0.0, keys(kr))
                splice!(trimmed, findfirst(ff -> ff.label == lbl, trimmed))
                vals[lbl] = uv(0.0, σ(kr, lbl))
            end
            if isempty(vals)
                return cat(append!(retained, removed))
            end
            push!(removed, uvs(vals))
        else
            return cat(retained)
        end
    end # while
    @assert false
    return cat(removed) # To maintain type
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
    forcezeros::Bool = true,
) where { T <: AbstractFloat }
    krs = _filterfit(unk, ffs, forcezeros)
    resid = Deferred() do
        sp = unk.label.spectrum
        props = copy(NeXLCore.properties(sp))
        props[:Name] = "Residual[$(props[:Name])]"
        residual = copy(sp.counts)
        for ff in ffs
            residual[ff.roi] -= (value(krs, ff.label) * ff.scale / unk.scale) * ff.charonly
        end
        return Spectrum(sp.energy, residual, props)
    end
    pb = Deferred() do 
        residual = copy(unk.label.spectrum.counts)
        for ff in ffs
            residual[ff.roi] -= (value(krs, ff.label) * ff.scale / unk.scale) * ff.charonly
        end
        res = Dict{ReferenceLabel,NTuple{3,T}}()
        for ff in ffs
            sr, sco = sum(residual[ff.roi]), sum(ff.charonly)
            res[ff.label] = (
                sr + sco*(value(krs, ff.label) * ff.scale / unk.scale), # All fitted counts
                sr, # Background only
                sco * ff.scale
            )
        end
        return res
    end
    return FilterFitResult(unk.label, krs, unk.roi, unk.data, resid, pb)
end

function fit_spectrum(
    ty::Type{FilteredUnknownW{T}},
    unk::Spectrum,
    filt::TopHatFilter,
    refs::AbstractVector{FilteredReference{T}},
    forcezeros::Bool = true,
) where { T <: AbstractFloat }
    bestRefs = selectBestReferences(refs)
    return filterfit(tophatfilter(ty, unk, filt, one(T) / T(dose(unk))), bestRefs, forcezeros)
end

function fit_spectra(
    ty::Type{FilteredUnknownW{T}},
    unks::AbstractVector{<:Spectrum},
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat }
    bestRefs = selectBestReferences(refs)
    return map(unks) do unk
        fu = tophatfilter(ty, unk, filt, one(T) / T(dose(unk)))
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

fit_spectra(
    unks::AbstractVector{Spectrum},
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true,
) where { T <: AbstractFloat } = # 
    fit_spectra(FilteredUnknownW{T}, unks, filt, refs, forcezeros)

fit_spectrum(
    unks::AbstractVector{Spectrum},
    filt::TopHatFilter{T},
    refs::AbstractVector{FilteredReference{T}},
    forcezeros = true) where { T<: AbstractFloat }= fit_spectra(unks, filt, refs, forcezeros)