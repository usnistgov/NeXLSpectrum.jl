using SparseArrays
using Polynomials
using LinearAlgebra
using NeXLUncertainties
using TimerOutputs

abstract type TopHatFilter end
struct VariableWidthFilter <: TopHatFilter end # Variable width, sparse filter
struct ConstantWidthFilter <: TopHatFilter end # Constant width, sparse filter

"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=1.0)::AbstractMatrix{Float64}

Build the default top-hat filter for the specified detector with the specified top and base parameters.
"""
buildfilter(det::Detector, a::Float64 = 1.0, b::Float64 = 1.0, full::Bool = false)::AbstractMatrix{Float64} =
    buildfilter(VariableWidthFilter, det, a, b, full)

"""
    buildfilter(::Type{<:TopHatFilter}, det::Detector, a::Float64=1.0, b::Float64=1.0)::AbstractMatrix{Float64}

Build a top-hat filter whose width varies with the detector's resolution as a function of X-ray energy for the
specified detector with the specified top and base parameters.
"""
function buildfilter(
    ty::Type{<:TopHatFilter},
    det::Detector,
    a::Float64 = 1.0,
    b::Float64 = 1.0
)::AbstractMatrix{Float64}
    resf(::Type{VariableWidthFilter}, ee, det) = resolution(ee, det)
    resf(::Type{ConstantWidthFilter}, ee, det) = resolution(energy(n"Mn K-L3"), det)
    intersect(a, b, c, d) = max(0.0, min(b, d) - max(a, c)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1) * (-0.5 / (mina - minb)) + intersect(mina, maxa, e0, e1) / (maxa - mina) +
        intersect(maxa, maxb, e0, e1) * (-0.5 / (maxb - maxa))
    cc = channelcount(det)
    ans = zeros(Float64, (cc, cc))
    for ch1 = 1:cc
        ee = 0.5 * (energy(ch1, det) + energy(ch1 + 1, det)) # midpoint of channel
        res = resf(ty, ee, det)
        ea = (ee - 0.5 * a * res, ee + 0.5 * a * res)
        eb = (ea[1] - 0.5 * b * res, ea[2] + 0.5 * b * res)
        chmin, chmax = channel(eb[1], det), channel(eb[2], det)
        if (chmin >= det.lld) && (chmax <= cc) # fit the full filter
            sum = 0.0 # Ensure the sum of the top-hat is precisely 0.0
            for ch2 = chmin:chmax-1
                fi = filtint(energy(ch2, det), energy(ch2 + 1, det), eb[1], ea[1], ea[2], eb[2])
                sum += fi
                ans[ch1, ch2] = fi
            end
            ans[ch1, chmax] = -sum # Ensure the sum is zero...
        end
    end
    return sparse(ans)
end

abstract type FilteredLabel <: Label end

struct ReferenceLabel <: FilteredLabel
    spec::Spectrum
    roi::UnitRange{Int}
end

Base.show(io::IO, refLab::ReferenceLabel) = print(io::IO, "$(refLab.spec[:Name])[$(refLab.roi)]")
Base.isequal(rl1::ReferenceLabel, rl2::ReferenceLabel) = isequal(rl1.roi, rl2.roi) && isequal(rl1.spec, rl2.spec)

struct UnknownLabel <: FilteredLabel
    spec::Spectrum
end

Base.show(io::IO, unk::UnknownLabel) = print(io, "Unknown[$(unk.spec[:Name])]")
Base.isequal(ul1::UnknownLabel, ul2::UnknownLabel) = isequal(ul1.spec, ul2.spec)

abstract type FilteredDatum end

"""
    FilteredReference

Represents the filtered reference spectrum over an ROI.
Carries the minimal data necessary to support filter-fitting a single
region-of-interest (continguous range of channles) and computing
useful output statistics.
"""
struct FilteredReference <: FilteredDatum
    identifier::ReferenceLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    back::Vector{Float64} # Background corrected intensity data over roi
end

Base.show(io::IO, fd::FilteredReference) = print(io, "Reference[$(fd.identifier)]")

"""
    FilteredUnknown

Represents the unknown in a filter fit.
"""
struct FilteredUnknown <: FilteredDatum
    identifier::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::AbstractMatrix{Float64} # Channel covariance
end

Base.show(io::IO, fd::FilteredUnknown) = print(io, "Unknown[$(fd.identifier)]")

function computeResidual(unk::FilteredUnknown, ffs::Array{FilteredReference}, kr::UncertainValues)
    res = copy(unk.data)
    for ff in filter(ff -> !ismissing(ff.back), ffs)
        ref = (NeXLUncertainties.value(ff.identifier, kr) * ff.scale / unk.scale) * ff.back
        for ch in ff.roi
            res[ch] -= ref[ch-ff.roi.start+1]
        end
    end
    return res
end

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::AbstractMatrix{Float64},
      back::Union{Missing,Vector{Float64}} = missing,
      scale = 1.0,
      tol::Float64 = 1.0e-4,
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum
with the specified filter.
"""
function filterImpl(
    spec::Spectrum,
    roi::UnitRange{Int},
    filter::AbstractMatrix{Float64},
    tol::Float64
)
    data = counts(spec, Float64) # Extract the spectrum data as Float64 to match the filter
    # Determine tangents to the two background end points
    tangents = map(st -> estimateBackground(data, st, 5, 2), (roi.start, roi.stop))
    # Replace the non-ROI channels with extensions of the tangent functions
    data[1:roi.start-1] = map(tangents[1], 1-roi.start:-1)
    data[roi.stop+1:end] = map(tangents[2], 1:length(data)-roi.stop)
    # Compute the filtered data
    filtered = filter * data
    maxval = maximum(filtered)
    roiff = findfirst(f -> abs(f) > tol * maxval, filtered):findlast(f -> abs(f) > tol * maxval, filtered)
    return ( roiff, data[roiff], filtered[roiff] )
end

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::AbstractMatrix{Float64},
      ashell::AtomicShell,
      scale = 1.0,
      tol = 1.0e-6
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.  Use a simple
edge-based background model.
"""
function Base.filter(
    spec::Spectrum,
    roi::UnitRange{Int},
    filt::AbstractMatrix{Float64},
    ashell::AtomicShell,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6
)::FilteredReference
    back = spec[roi] - modelBackground(spec, roi, ashell)
    f = filterImpl(spec, roi, filt, tol)
    return FilteredReference(ReferenceLabel(spec, roi), scale, roi, f..., back)
end

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::AbstractMatrix{Float64},
      scale = 1.0,
      tol = 1.0e-6
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter. Use a naive
linear background model.
"""
function Base.filter(
    spec::Spectrum,
    roi::UnitRange{Int},
    filt::AbstractMatrix{Float64},
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6
)::FilteredReference
    back = spec[roi] - modelBackground(spec, roi)
    f = filterImpl(spec, roi, filt, tol)
    return FilteredReference(ReferenceLabel(spec, roi), scale, roi, f..., back)
end

"""
    filter(spec::Spectrum, filter::AbstractMatrix{Float64}, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredDatum

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
function Base.filter(spec::Spectrum, filt::AbstractMatrix{Float64}, scale::Float64 = 1.0)::FilteredUnknown
    data = counts(spec, Float64) # Extract the spectrum data
    # Compute the filtered data
    filtered = filt * data
    roi = 1:findlast(f -> f ≠ 0.0, filtered)
    # max(d,1.0) is necessary to ensure the variances are positive
    diag = sparse(roi,roi,map(d -> max(d, 1.0), data[roi]))
    # diag = Diagonal(map(d -> max(d, 1.0), data[roi]))
    covar = filt[roi,roi] * diag * transpose(filt[roi,roi])
    @assert covar isa AbstractSparseMatrix "The covariance matrix should be a sparse matrix in filter(...)::FilterUnknown."
    checkcovariance!(covar)
    return FilteredUnknown(UnknownLabel(spec), scale, roi, roi, data[roi], filtered[roi], covar)
end

"""
    extract(fd::FilteredReference, roi::UnitRange{Int})

Extract the filtered data representing the specified range.  <code>roi</code> must fully encompass the filtered
data in <code>fd</code>.
"""
function extract(fd::FilteredReference, roi::UnitRange{Int})
    @assert fd.ffroi.start>=roi.start "$(fd.ffroi.start) < $(roi.start) in $(fd)"
    @assert fd.ffroi.stop<=roi.stop "$(fd.ffroi.stop) > $(roi.stop) in $(fd)"
    data = zeros(Float64, length(roi))
    nz = fd.ffroi.start-roi.start+1:fd.ffroi.stop-roi.start+1
    data[nz] = fd.filtered
    return data
end

"""
    extract(fd::FilteredUnknown, roi::UnitRange{Int})

Extract the filtered data representing the specified range.  <code>roi</code> must fully encompass the filtered
data in <code>fd</code>.
"""
function extract(fd::FilteredUnknown, roi::UnitRange{Int})
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return copy(fd.filtered[nz])
end

"""
    covariance(fd::FilteredUnknown, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.  <code>roi</code> must
fully encompass the filtered edata in <code>fd</code>.
"""
function covariance(fd::FilteredUnknown, roi::UnitRange{Int})
    # Unknown case
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return Matrix(fd.covariance[nz,nz])
end

"""
Generalized least squares (my implementation)
"""
function fitcontiguousg(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * glssvd(extract(unk, chs), A, covariance(unk, chs), xlbls)
end

"""
Generalized least squares (pseudo-inverse)
"""
function fitcontiguousp(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * glspinv(extract(unk, chs), A, covariance(unk, chs), xlbls)
end

"""
Weighted least squares
"""
function fitcontiguousw(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * wlssvd(extract(unk, chs), A, diag(covariance(unk, chs)), xlbls)
end

"""
Ordinary least squares
"""
function fitcontiguouso(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * olssvd(extract(unk, chs), A, 1.0, xlbls)
end

struct FilterFitResult
    label::UnknownLabel
    kratios::UncertainValues
    roi::UnitRange{Int}
    raw::Vector{Float64}
    residual::Vector{Float64}
end

tabulate(ffrs::Array{FilterFitResult}, withUnc=false) =
    NeXLUncertainties.tabulate([r.kratios for r in ffrs],withUnc)

Base.show(io::IO, ffr::FilterFitResult) = print(io, "$(ffr.label)")

NeXLUncertainties.value(label::ReferenceLabel, ffr::FilterFitResult) = value(label, ffr.kratios)

NeXLUncertainties.σ(label::ReferenceLabel, ffr::FilterFitResult) = σ(label, ffr.kratios)

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
    filterfit(unk::FilteredUnknown, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknown, ffs::Array{FilteredReference}, alg = fitcontiguousp)::FilterFitResult
    fitrois = ascontiguous(map(fd->fd.ffroi, ffs))
    kr = cat(map(fr->alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs), fr), fitrois))
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, computeResidual(unk, ffs, kr))
end

"""
    filteredresidual(fit::UncertainValues, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
function filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}
    scaled(ff) = (value(ff.identifier, fit) * (ff.scale / unk.scale)) * extract(ff, unk.ffroi)
    return unk.filtered - mapreduce(scaled, +, ffs)
end
