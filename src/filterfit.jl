using SparseArrays
using Polynomials
using LinearAlgebra
using NeXLUncertainties
using TimerOutputs

"""
The TopHatFilter struct represents a zero-sum symmetric second-derivative-like filter that when applied to
spectral data has the property of suppressing constant and slowly varying signals (like the continuum) while
retaining a linear signal for faster changing signals like the characteristic peaks.

See
  * F. H. Schamber Proc Symposium of "X-ray Fluorscence Analysis on Environmental Samples" Chapel Hill 1976 T Dzubay Ed.
  * P. Statham Anal Chem 49 no 14 Dec 1977

The TopHatFilter struct optimizes the memory and CPU use when applying top-hat filters to spectrum data.

The easiest way to implement a top-hat filter is as matrix F.  The rows represent the filters.  The product of the
filter and the data vector is the filtered spectrum.  The product of the filter times a diagnonal matrix constructed
from the data times the transpose of the filter is the covariance of the filtered data.  The diagonal matrix constructed
from the spectrum data is the covariance matrix associated with the spectrum data because the channels in the spectrum
data are independent (thus the matrix is diagnonal) and the magnitude equals the counts in each channels because the
spectrum data is nominally Poissonian and in the large number limit, the variance of a Poissonian random variable is
the number itself (σ=sqrt(N) => Var = N)

The filter matrix however is sparse and we can use the SparseMatrix type to implement it with better memory and CPU
performance.  But we can do better. The filter matrix mostly zeros except for a small band near the diagonal.  This
band is continuous and can readily be implemented as an array.  Furthermore, the product F⋅d and F⋅D⋅Fᵀ are readily
implemented as element-by-element multiplies and dot products.  Thus storing the filter as offsets and row filters is
efficient in both memory and CPU.
"""
struct TopHatFilter
    offsets::Vector{Int}  # offset to start of filter in row
    filters::Vector{Vector{Float64}} # Filter data
    weights::AbstractVector{Float64} # Correlation compensation weights

    function TopHatFilter(filt::AbstractMatrix{Float64}, wgts::AbstractVector{Float64})
        # Extract the contiguous non-zero elements out
        dim = size(filt)[1]
        offsets = zeros(dim)
        filts = fill(Vector{Float64}(), dim)
        for (r, row) in enumerate(eachrow(filt))
            start = findfirst(i -> i ≠ 0.0, row)
            if !isnothing(start)
                stop = findlast(i -> i ≠ 0.0, row)
                offsets[r] = start
                filts[r] = [ row[start:stop]... ]
            end
        end
        return new(offsets, filts, wgts)
    end
end

abstract type FittingFilterType end

struct VariableWidthFilter <: FittingFilterType end # Variable width, sparse filter
struct ConstantWidthFilter <: FittingFilterType end # Constant width, sparse filter


"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=1.0)::TopHatFilter

Build the default top-hat filter for the specified detector with the specified top and base parameters.
"""
buildfilter(det::Detector, a::Float64 = 1.0, b::Float64 = 1.0, full::Bool = false)::TopHatFilter =
    buildfilter(VariableWidthFilter, det, a, b, full)

"""
    buildfilter(::Type{<:FittingFilterType}, det::Detector, a::Float64=1.0, b::Float64=1.0)::TopHatFilter

Build a top-hat filter whose width varies with the detector's resolution as a function of X-ray energy for the
specified detector with the specified top and base parameters.
"""
function buildfilter(
    ty::Type{<:FittingFilterType},
    det::Detector,
    a::Float64 = 1.0, # Top
    b::Float64 = 1.0  # Base
)::TopHatFilter
    resf(::Type{VariableWidthFilter}, ee, det) = resolution(ee, det) # FWHM at ee
    resf(::Type{ConstantWidthFilter}, ee, det) = resolution(energy(n"Mn K-L3"), det) # FWHM at Mn Ka
    intersect(a, b, c, d) = max(0.0, min(b, d) - max(a, c)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1) * (-0.5 / (mina - minb)) + intersect(mina, maxa, e0, e1) / (maxa - mina) +
        intersect(maxa, maxb, e0, e1) * (-0.5 / (maxb - maxa))
    M(minb, mina) = (channel(mina, det) - channel(minb, det))/2 # base length in channels
    N(mina, maxa) = (channel(mina, det) - channel(maxa, det))/2 # Top length in channels
    cc = channelcount(det)
    filt, wgts = zeros(Float64, (cc, cc)), zeros(Float64, cc)
    for ch1 in eachindex(wgts)
        ee = 0.5 * (energy(ch1, det) + energy(ch1 + 1, det)) # midpoint of channel
        res = resf(ty, ee, det)
        ea = ( ee - 0.5 * a * res, ee + 0.5 * a * res )
        eb = ( ea[1] - 0.5 * b * res, ea[2] + 0.5 * b * res )
        chmin, chmax = channel(eb[1], det), channel(eb[2], det)
        if true
            if (chmin>=det.lld) && (chmax<=cc)
                for ch2 in chmin:chmax
                    filt[ch1, ch2] = filtint(energy(ch2, det), energy(ch2 + 1, det), eb[1], ea[1], ea[2], eb[2])
                end
            end
        else
            for ch2 in max(chmin, det.lld):min(chmax, cc)
                filt[ch1, ch2] = filtint(energy(ch2, det), energy(ch2 + 1, det), eb[1], ea[1], ea[2], eb[2])
            end
        end
        # Ensure the sum is zero...
        if chmin<det.lld
            filt[ch1,det.lld] = -sum(filt[ch1, det.lld+1:chmax])
        else
            filt[ch1, min(cc,chmax)] = -sum(filt[ch1,chmin:min(cc,chmax)-1])
        end
        m, n = M(eb[1], ea[1]), N(ea...)
        wgts[ch1] = ((2.0*m+1.0)*(2.0*n))/(2.0*m+1.0+2.0*n) # Schamber's formula (see Statham 1977, Anal Chem 49, 14 )
    end
    return TopHatFilter(filt, wgts)
end

abstract type FilteredLabel <: Label end

abstract type ReferenceLabel <: FilteredLabel end

struct SpectrumLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
end


Base.show(io::IO, refLab::SpectrumLabel) = print(io::IO, "$(refLab.spec[:Name])[$(refLab.roi)]")
Base.isequal(rl1::SpectrumLabel, rl2::SpectrumLabel) = isequal(rl1.roi, rl2.roi) && isequal(rl1.spec, rl2.spec)

struct CharXRayLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
    xrays::Vector{CharXRay}
end

Base.show(io::IO, refLab::CharXRayLabel) = print(io::IO, "$(name(refLab.xrays))")
Base.isequal(rl1::CharXRayLabel, rl2::CharXRayLabel) =
    isequal(rl1.roi, rl2.roi) && isequal(rl1.xrays, rl2.xrays) && isequal(rl1.spec, rl2.spec)

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
    covscale::Float64  # Schamber's variance scale factor
end

Base.show(io::IO, fd::FilteredReference) = print(io, "Reference[$(fd.identifier)]")

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::TopHatFilter,
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
    filter::TopHatFilter,
    tol::Float64
)
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
    data = counts(spec, Float64) # Extract the spectrum data as Float64 to match the filter
    # Determine tangents to the two background end points
    tangents = map(st -> estimatebackground(data, st, 5, 2), (roi.start, roi.stop))
    # Replace the non-ROI channels with extensions of the tangent functions
    data[1:roi.start-1] = map(tangents[1], 1-roi.start:-1)
    data[roi.stop+1:end] = map(tangents[2], 1:length(data)-roi.stop)
    # Compute the filtered data
    filtered = apply(filter, data)
    maxval = maximum(filtered)
    roiff = findfirst(f -> abs(f) > tol * maxval, filtered):findlast(f -> abs(f) > tol * maxval, filtered)
    return ( roiff, data[roiff], filtered[roiff] )
end

function Base.filter(
    spec::Spectrum,
    det::Detector,
    cxrs::AbstractVector{CharXRay},
    filter::TopHatFilter,
    scale::Float64=1.0,
    ampl::Float64=1.0e-4,
    tol::Float64=1.0e-6
)::Vector{FilteredReference}
    res=[]
    for (lbl, roi) in labeledextents(cxrs, det, ampl)
        back = spec[roi] - modelBackground(spec, roi, inner(brightest(lbl)))
        f = filterImpl(spec, roi, filter, tol)
        push!(res, FilteredReference(CharXRayLabel(spec, roi, lbl), scale, roi, f..., back, filter.weights[(roi.start+roi.stop)÷2]))
    end
    return res
end

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::TopHatFilter,
      ashell::AtomicSubShell,
      scale = 1.0,
      tol = 1.0e-6
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.  Use a simple
edge-based background model.
"""
function Base.filter(
    spec::Spectrum,
    roi::UnitRange{Int},
    filter::TopHatFilter,
    ashell::AtomicSubShell,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6
)::FilteredReference
    back = spec[roi] - modelBackground(spec, roi, ashell)
    f = filterImpl(spec, roi, filter, tol)
    return FilteredReference(SpectrumLabel(spec, roi), scale, roi, f..., back, filter.weights[(roi.start+roi.stop)÷2])
end

"""
    filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::TopHatFilter,
      scale = 1.0,
      tol = 1.0e-6
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter. Use a naive
linear background model.
"""
function Base.filter(
    spec::Spectrum,
    roi::UnitRange{Int},
    filter::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6
)::FilteredReference
    back = spec[roi] - modelBackground(spec, roi)
    f = filterImpl(spec, roi, filter, tol)
    return FilteredReference(SpectrumLabel(spec, roi), scale, roi, f..., back, filter.weights[(roi.start+roi.stop)÷2])
end


abstract type FilteredUnknown <: FilteredDatum end


"""
    FilteredUnknownG

Represents the unknown in a filter fit using the full generalized fitting model.
"""
struct FilteredUnknownG <: FilteredUnknown
    identifier::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::AbstractMatrix{Float64} # Channel covariance
end

Base.show(io::IO, fd::FilteredUnknown) = print(io, "Unknown[$(fd.identifier)]")

"""
    filter(spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredDatum

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
function Base.filter(::Type{FilteredUnknownG}, spec::Spectrum, filter::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownG
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
    function covariance(roi, data, filter)
        off(r, o) = r.start-o+1:r.stop-o+1
        rs, cs, vs = Vector{Int}(), Vector{Int}(), Vector{Float64}()
        res = zeros( ( length(roi), length(roi) ) )
        bdata = map(d -> max(d, 1.0), data)
        for r in roi
            rr, rf, ro = range(r), filter.filters[r], filter.offsets[r]
            for c in roi
                ii = intersect(rr, range(c))
                if length(ii) > 0
                    cf, co = filter.filters[c], filter.offsets[c]
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
    data = counts(spec, Float64) # Extract the spectrum data
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
    covariance(fd::FilteredUnknown, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.  <code>roi</code> must
fully encompass the filtered edata in <code>fd</code>.
"""
function covariance(fd::FilteredUnknownG, roi::UnitRange{Int})
    # Unknown case
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return Matrix(fd.covariance[nz,nz])
end

"""
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.
"""
struct FilteredUnknownW <: FilteredUnknown
    identifier::UnknownLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::Vector{Float64} # Channel covariance
end

"""
    filter(spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknown

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
Base.filter(spec::Spectrum, filt::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknown =
    filter(FilteredUnknownG,spec,filt,scale)


"""
    filter(::Type{FilteredUnknownW}, spec::Spectrum, filter::TopHatFilter, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredUnknownW

For filtering the unknown spectrum. Process the full Spectrum with the specified filter for use with the weighted
least squares model.
"""
function Base.filter(::Type{FilteredUnknownW}, spec::Spectrum, filter::TopHatFilter, scale::Float64 = 1.0)::FilteredUnknownW
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
    covariance(filter, data) =
        [ dot(filter.filters[r] .* data[range(r)], filter.filters[r]) for r in eachindex(data) ]
    data = counts(spec, Float64) # Extract the spectrum data
    # Compute the filtered data
    filtered = apply(filter, data)
    roi = 1:findlast(f -> f ≠ 0.0, filtered)
    # max(d,1.0) is necessary to ensure the variances are positive
    covar = covariance(filter, map(d -> max(d, 1.0), data))
    return FilteredUnknownW(UnknownLabel(spec), scale, roi, roi, data[roi], filtered[roi], covar)
end


"""
    covariance(fd::FilteredUnknown, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance diagnonal elements over the specified range of channels.
<code>roi</code> must fully encompass the filtered edata in <code>fd</code>.
"""
function covariance(fd::FilteredUnknownW, roi::UnitRange{Int})
    # Unknown case
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return fd.covariance[nz]
end

"""
    computeResidual(unk::FilteredUnknown, ffs::Array{FilteredReference}, kr::UncertainValues)

Computes the residual spectrum for the specified unknown.
"""

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
Generalized least squares (my implementation)
"""
function fitcontiguousg(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * glssvd(extract(unk, chs), x, covariance(unk, chs), lbls)
end

"""
Generalized least squares (pseudo-inverse)
"""
function fitcontiguousp(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * glspinv(extract(unk, chs), x, covariance(unk, chs), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * glspinv(extract(unk, chs), x, Diagonal(diag(covariance(unk, chs))), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw2(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    return scale * wlssvd(extract(unk, chs), x, diag(covariance(unk, chs)), lbls)
end

"""
Weighted least squares for FilteredUnknownW
"""
function fitcontiguousww(unk::FilteredUnknownW, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    wgts = [ ff.covscale for ff in ffs ]
    # return scale * wlssvd(extract(unk, chs), A, covariance(unk, chs), xlbls)
    return scale * wlspinv(extract(unk, chs), x, covariance(unk, chs), wgts, lbls)
end

"""
Ordinary least squares for either FilteredUnknown[G|W]
"""
function fitcontiguouso(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    # Build labels and scale
    lbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])
    # return scale * olssvd(extract(unk, chs), A, 1.0, xlbls)
    return scale * olspinv(extract(unk, chs), x, 1.0, lbls)
end


"""
    FilterFitResult

Represents the result of fitting either FilteredUnknownG or FilteredUnknownW to a FilteredUnknown.

Struct elements

    label::UnknownLabel  # Identifies the unknown
    kratios::UncertainValues # Labeled with ReferenceLabel objects
    roi::UnitRange{Int} # Range of channels fit
    raw::Vector{Float64} # Raw spectrum data
    residual::Vector{Float64} # Residual spectrum
"""
struct FilterFitResult
    label::UnknownLabel
    kratios::UncertainValues
    roi::UnitRange{Int}
    raw::Vector{Float64}
    residual::Vector{Float64}
end

"""
    kratios(ffr::FilterFitResult)

The k-ratios as a UncertainValues object
"""
kratios(ffr::FilterFitResult)::UncertainValues = ffr.kratios
unknown(ffr::FilterFitResult)::UnknownLabel = ffr.label
residual(ffr::FilterFitResult)::Vector{Float64} = ffr.residual

tabulate(ffrs::Array{FilterFitResult}, withUnc=false) =
    NeXLUncertainties.tabulate([r.kratios for r in ffrs],withUnc)

function details(io::IO, ffr::FilterFitResult)
    println(io, "  Unknown:   $(ffr.label)")
    for l in labels(ffr)
        println(io, "   $(l): $(ffr[l])")
    end
end

details(ffr::FilterFitResult) = details(stdout,ffr)

Base.show(io::IO, ffr::FilterFitResult) = print(io, "$(ffr.label)")

NeXLUncertainties.value(label::ReferenceLabel, ffr::FilterFitResult) = value(label, ffr.kratios)
NeXLUncertainties.σ(label::ReferenceLabel, ffr::FilterFitResult) = σ(label, ffr.kratios)
NeXLUncertainties.getindex(ffr::FilterFitResult, label::ReferenceLabel) = getindex(ffr.kratios, label)
NeXLUncertainties.labels(ffr::FilterFitResult) = labels(ffr.kratios)

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
    filterfit(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg = fitcontiguousp)::FilterFitResult
    fitrois = ascontiguous(map(fd->fd.ffroi, ffs))
    kr = cat(map(fr->alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs), fr), fitrois))
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, computeResidual(unk, ffs, kr))
end

"""
    filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg = fitcontiguousww)::FilterFitResult
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
