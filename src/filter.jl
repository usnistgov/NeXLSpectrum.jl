using SparseArrays
using Polynomials
using LinearAlgebra

"""
    FittingFilterType

Represents different fitting filter models.  A fitting filter is a linear operator that is 1) symmetric about the
center, and 2) the sum of the elements must be zero.  The standard filter is a top-hat shape, although Gaussian,
triangular, Savitsky-Golay and other shapes are also possible.
"""
abstract type FittingFilterType end

"""
    VariableWidthFilter

A top-hat filter that varies in width with the FWHM of the detector.
"""
struct VariableWidthFilter <: FittingFilterType end # Variable width filter

"""
    ConstantWidthFilter

A top-hat filter that has constant width determined by FWHM at Mn Kα for all channels.
"""
struct ConstantWidthFilter <: FittingFilterType end # Constant width filter

"""
    GaussianFilter

A Gaussian-shaped filter that varies in width with the FWHM of the detector.  The Gaussian is offset to ensure the
sum of the filter elements is zero.
"""
struct GaussianFilter <: FittingFilterType end # A variable width Gaussian filter

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
    filttype::Type{<:FittingFilterType}
    detector::Detector
    offsets::Vector{Int}  # offset to start of filter in row
    filters::Vector{Vector{Float64}} # Filter data
    weights::AbstractVector{Float64} # Correlation compensation weights

    function TopHatFilter(
        ty::Type{<:FittingFilterType},
        det::Detector,
        filt::AbstractMatrix{Float64},
        wgts::AbstractVector{Float64},
    )
        # Extract the contiguous non-zero elements out
        dim = size(filt)[1]
        offsets = zeros(dim)
        filts = fill(Vector{Float64}(), dim)
        for (r, row) in enumerate(eachrow(filt))
            start = findfirst(i -> i ≠ 0.0, row)
            if !isnothing(start)
                stop = findlast(i -> i ≠ 0.0, row)
                offsets[r] = start
                filts[r] = [row[start:stop]...]
            end
        end
        return new(ty, det, offsets, filts, wgts)
    end
end

function filterdata(filt::TopHatFilter, row::Int)::Vector{Float64}
    res = zeros(Float64, length(filt.filters))
    res[filt.offsets[row]:filt.offsets[row]+length(filt.filters[row])-1] = filt.filters[row]
    return res
end

Base.size(filt::TopHatFilter) = ( length(filt.filters), length(filt.filters) )


Base.show(io::IO, thf::TopHatFilter) = print(io, "$(thf.filttype)[$(thf.detector)]")

"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=2.0)::TopHatFilter

Build the default top-hat filter for the specified detector with the specified top and base parameters.
"""
buildfilter(det::Detector, a::Float64 = 1.0, b::Float64 = 2.0)::TopHatFilter =
    buildfilter(VariableWidthFilter, det, a, b)

"""
    buildfilter(::Type{<:FittingFilterType}, det::Detector, a::Float64=1.0, b::Float64=1.0)::TopHatFilter

Build a top-hat-style filter for the specified detector with the specified top and base parameters. The
VariableWidthFilter and ConstantWidthFilter types are currently supported.
"""
function buildfilter(
    ty::Type{<:FittingFilterType},
    det::Detector,
    a::Float64 = 1.0, # Top
    b::Float64 = 2.0,  # Base
)::TopHatFilter
    resf(::Type{VariableWidthFilter}, center, det) = resolution(center, det) # FWHM at center
    resf(::Type{ConstantWidthFilter}, center, det) = resolution(energy(n"Mn K-L3"), det) # FWHM at Mn Ka
    intersect(aa, bb, cc, dd) = max(0.0, min(bb, dd) - max(aa, cc)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1) * (-0.5 / (mina - minb)) + intersect(mina, maxa, e0, e1) / (maxa - mina) +
        intersect(maxa, maxb, e0, e1) * (-0.5 / (maxb - maxa))
    M(astop, astart) = (channel(astop, det) - channel(astart, det) - 1.0) / 2.0   # base length in channels
    N(bstart, astart) = channel(astart, det) - channel(bstart, det) # Top length in channels
    cc = channelcount(det)
    filt, wgts = zeros(Float64, (cc, cc)), zeros(Float64, cc)
    for ch1 in eachindex(wgts)
        center = 0.5 * (energy(ch1, det) + energy(ch1 + 1, det)) # midpoint of channel
        res = resf(ty, center, det)
        ea = (center - 0.5 * a * res, center + 0.5 * a * res)
        eb = (ea[1] - 0.5 * b * res, ea[2] + 0.5 * b * res)
        chmin, chmax = channel(eb[1], det), channel(eb[2], det)
        if (chmin >= 1) && (chmax <= cc)
            for i = 0:(chmax-chmin)÷2
                filt[ch1, chmax-i] = (filt[ch1, i+chmin] = filtint(
                    energy(i + chmin, det),
                    energy(i + chmin + 1, det),
                    eb[1],
                    ea[1],
                    ea[2],
                    eb[2],
                ))
            end
            filt[ch1, chmax] = (filt[ch1, chmin] -= sum(filt[ch1, chmin:min(cc, chmax)]) / 2.0)
            @assert abs(sum(filt[ch1, :])) < 1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i -> filt[ch1, i] == filt[ch1, chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - O"
        end
        m, n = M(ea[2], ea[1]), N(eb[1], ea[1])
        @assert m > 0.0 "m=$(m)!!!!"
        @assert n > 0.0 "n=$(n)!!!!"
        wgts[ch1] = 3.0 * ((2.0 * m + 1.0) * (2.0 * n)) / (2.0 * m + 1.0 + 2.0 * n) # Schamber's formula (see Statham 1977, Anal Chem 49, 14 )
    end
    @assert all(r -> abs(sum(r)) < 1.0e-12, eachrow(filt)) "A filter does not sum to zero."
    return TopHatFilter(ty, det, filt, wgts)
end

"""
    buildfilter(::Type{GaussianFilter}, det::Detector, a::Float64=1.0, b::Float64=5.0)::TopHatFilter

Build a top-hat filter with Gaussian shape whose width varies with the detector's resolution as a function of X-ray
energy for the specified detector with the specified top and base parameters. The <code>a</code> parameter corresponds
to the filter width relative to the detector resolution expressed as Gaussian width.  So <code>a=1</code> is a filter
whose width equals the detector resolution at each energy.  The <code>b</code> parameter is the extent of the
filter in Gaussian widths.  The default <code>a=1, b=5</code> corresponds to a  filter that has the same resolution
as the detector and an extent of 2.5 Gaussian widths above and below the center channel.
"""
function buildfilter(
    ty::Type{GaussianFilter},
    det::Detector,
    a::Float64 = 1.0,
    b::Float64 = 6.0, # Width
)::TopHatFilter
    filtint(center, ee, gw) = exp(-0.5 * ((ee - center) / gw)^2)
    cc = channelcount(det)
    filt, wgts = zeros(Float64, (cc, cc)), zeros(Float64, cc)
    for ch1 in eachindex(wgts)
        center = energy(ch1, det) # midpoint of channel
        res = a * gaussianwidth(resolution(center, det))
        ex = (center - 0.5 * b * res, center + 0.5 * b * res)
        chmin, chmax = channel(ex[1], det), channel(ex[2], det)
        if (chmin >= 1) && (chmax <= cc)
            for i = 0:(chmax-chmin)÷2 # Ensure that it is symmetric
                filt[ch1, chmax-i] = (filt[ch1, chmin+i] = filtint(center, energy(chmin + i, det), res))
            end
            # Offset the Gaussian to ensure the sum is zero.
            filt[ch1, chmin:chmax] .-= sum(filt[ch1, chmin:chmax]) / length(chmin:chmax)
            @assert abs(sum(filt[ch1, :])) < 1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i -> filt[ch1, i] == filt[ch1, chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - G"
        end
        wgts[ch1] = 1.0
    end
    return TopHatFilter(ty, det, filt, wgts)
end

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
    charonly::Vector{Float64} # Background corrected intensity data over roi
    covscale::Float64  # Schamber's variance scale factor
end

Base.show(io::IO, fd::FilteredReference) = print(io, "Reference[$(fd.identifier)]")

"""
    _filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      filter::TopHatFilter,
      tol::Float64
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.
"""
function _filter(spec::Spectrum, roi::UnitRange{Int}, filter::TopHatFilter, tol::Float64)
    range(i) = filter.offsets[i]:filter.offsets[i]+length(filter.filters[i])-1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [dot(filter.filters[i], view(data, range(i))) for i in eachindex(data)]
    # Extract the spectrum data as Float64 to match the filter
    data = counts(spec, Float64, true)
    # Determine tangents to the two background end points
    tangents = map(st -> estimatebackground(data, st, 5, 2), (roi.start, roi.stop))
    # Replace the non-ROI channels with extensions of the tangent functions
    data[1:roi.start-1] = map(tangents[1], 1-roi.start:-1)
    data[roi.stop+1:end] = map(tangents[2], 1:length(data)-roi.stop)
    # Compute the filtered data
    filtered = apply(filter, data)
    maxval = maximum(filtered)
    ff = findfirst(f -> abs(f) > tol * maxval, filtered)
    if isnothing(ff)
        error("There doesn't appear to be any data in $(spec[:Name]))")
    end
    roiff = ff:findlast(f -> abs(f) > tol * maxval, filtered)
    return (roiff, data[roiff], filtered[roiff])
end

Base.convert(::AbstractVector{SpectrumFeature}, x::AbstractVector{CharXRay})::AbstractVector{SpectrumFeature} = x

charFeature(elm::Element, tr::Tuple{Vararg{Transition}}; minweight = 1.0e-3, maxE = 1.0e6) =
    convert(AbstractVector{SpectrumFeature}, characteristic(elm, tr, minweight, maxE))

escapeFeature(elm::Element, trs::Tuple{Vararg{Transition}}; minweight = 0.1, maxE = 1.0e6, escape = n"Si K-L3") =
    convert(
        AbstractVector{SpectrumFeature},
        map(tr -> EscapeArtifact(tr, escape), characteristic(elm, trs, minweight, maxE)),
    )

comptonFeature(elm::Element, trs::Tuple{Vararg{Transition}}; minweight = 1.0e-3) = convert(
    AbstractVector{SpectrumFeature},
    map(tr -> ComptonArtifact(tr, escape), characteristic(elm, trs, minweight, maxE)),
)

function Base.filter(
    escLabel::EscapeLabel,
    filter::TopHatFilter,
    scale::Float64 = 1.0;
    ampl::Float64 = 1.0e-4,
)::FilteredReference
    spec, roi = esc.spec, esc.roi
    charonly = spec[roi] - modelBackground(spec, roi)
    f = _filter(spec, roi, filter, tol)
    return FilteredReference(
        escLabel,
        scale,
        roi,
        f...,
        charonly,
        filter.weights[(roi.start+roi.stop)÷2],
    )
end


"""
    function labeledextents(
        elm::Element,  # All CharXRay for this element
        det::Detector,
        ampl::Float64,
        maxE::Float64=energy(det.channelcount+1, det) # full detector range
    )::Vector{Tuple{Vector{CharXRay},UnitRange{Int}}}

Creates a vector containing pairs containing a vector of CharXRay and an interval. The interval represents a
contiguous interval over which all the X-rays in the interval are sufficiently close in energy that they will
interfere with each other on the specified detector.
"""
labeledextents(elm::Element, det::Detector, ampl::Float64, maxE::Float64=1.0e6) =
    labeledextents(visible(characteristic(elm, alltransitions, ampl, maxE), det), det, ampl)

"""
    charXRayLabels(#
      spec::Spectrum, #
      elm::Element, #
      allElms::Vector{Element}, #
      det::Detector, #
      ampl::Float64, #
      maxE::Float64=1.0e6)::Vector{CharXRayLabel}

Creates a vector CharXRayLabel objects associated with 'elm' for a spectrum containing the elements
'allElms' assuming that it was collected on 'det'.  ROIs in which other elements from 'allElms'
interfere with 'elm' will not be included.
"""
function charXRayLabels(#
    spec::Spectrum, #
    elm::Element, #
    allElms::Vector{Element}, #
    det::Detector, #
    ampl::Float64, #
    maxE::Float64=1.0e6)::Vector{CharXRayLabel}
    @assert elm in allElms "$elm must be in $allElms."
    intersects(urs, roi) = !isnothing(findfirst(ur->length(intersect(ur,roi))>0, urs))
    # Find all the ROIs associated with other elements
    urs = collect(Iterators.flatten(extents(ae, det, ampl) for ae in filter(a->a ≠ elm, allElms)))
    # Find elm's ROIs that don't intersect another element's ROI
    lxs = filter(lx->!intersects(urs, lx[2]), labeledextents(elm, det, ampl, maxE))
    return [ CharXRayLabel(spec,roi,xrays) for (xrays, roi) in lxs ]
end

"""
    Base.filter(
        charLabel::CharXRayLabel,
        filter::TopHatFilter,
        scale::Float64 = 1.0,
        tol::Float64 = 1.0e-6,
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.  Use a simple
edge-based background model.
"""
function Base.filter(
    charLabel::CharXRayLabel,
    filter::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::FilteredReference
    ashellof(xrays) =
        inner(xrays[argmax(jumpratio.(inner.(xrays)))])
    spec, roi, ashell = charLabel.spec, charLabel.roi, ashellof(charLabel.xrays)
    return FilteredReference(
        charLabel,
        scale,
        roi,
        _filter(spec, roi, filter, tol)...,
        spec[roi] - modelBackground(spec, roi, ashell),
        filter.weights[(roi.start+roi.stop)÷2],
    )
end

Base.filter(
    labels::Vector{<:ReferenceLabel},
    filt::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::Vector{FilteredReference} =
    map(lbl->filter(lbl,filt,scale,tol), labels)

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
    reflabel::ReferenceLabel,
    filter::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::FilteredReference
    spec, roi = reflabel.spec, reflabel.roi
    return FilteredReference(
        reflabel,
        scale,
        roi,
        _filter(spec, roi, filter, tol)...,
        spec[roi] - modelBackground(spec, roi),
        filter.weights[(roi.start+roi.stop)÷2],
    )
end

"""
    FilteredUnknown

A FilteredDatum representing the unknown.  There are two types of FilteredUnknown - FilteredUnknownG and
FilteredUnknownW for "generalized" and "weighted" model unknowns.  The distinction is desirable since the
full covariance calculation for the generalized model is so expensive relative to the weighted approximation.
"""
abstract type FilteredUnknown <: FilteredDatum end

Base.show(io::IO, fd::FilteredUnknown) = print(io, fd.identifier)


"""
    extract(fd::FilteredReference, roi::UnitRange{Int})

Extract the filtered data representing the specified range.  <code>roi</code> must fully encompass the filtered
data in <code>fd</code>.
"""
function extract(fd::FilteredReference, roi::UnitRange{Int})
    @assert fd.ffroi.start >= roi.start "$(fd.ffroi.start) < $(roi.start) in $(fd)"
    @assert fd.ffroi.stop <= roi.stop "$(fd.ffroi.stop) > $(roi.stop) in $(fd)"
    data = zeros(Float64, length(roi))
    nz = fd.ffroi.start-roi.start+1:fd.ffroi.stop-roi.start+1
    data[nz] = fd.filtered
    return data
end

"""
    extract(fd::FilteredUnknown, roi::UnitRange{Int})::AbstractVector{Float64}

Extract the filtered data representing the specified range.  <code>roi</code> must be fully contained within the
filtered data in <code>fd</code>.
"""
extract(fd::FilteredUnknown, roi::UnitRange{Int})::AbstractVector{Float64} = fd.filtered[roi]

_buildlabels(ffs::Array{FilteredReference}) = collect(ff.identifier for ff in ffs)
_buildscale(unk::FilteredUnknown, ffs::Array{FilteredReference}) =
    Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])

# Internal: Computes the residual spectrum based on the fit k-ratios
function _computeResidual(unk::FilteredUnknown, ffs::Array{FilteredReference}, kr::UncertainValues)
    res = copy(unk.data)
    for ff in ffs
        res[ff.roi] -= (value(ff.identifier, kr) * ff.scale / unk.scale) * ff.charonly
    end
    return res
end

# Internal: Computes the peak and background count based on the fit k-ratios
function _computecounts( #
    unk::FilteredUnknown, #
    ffs::Array{FilteredReference}, #
    kr::UncertainValues, #
)::Dict{<:ReferenceLabel,NTuple{2,Float64}}
    res = Dict{ReferenceLabel,NTuple{2,Float64}}()
    for ff in ffs
        su = sum(unk.data[ff.roi])
        res[ff.identifier] = (su, su - sum((value(ff.identifier, kr) * ff.scale / unk.scale) * ff.charonly))
    end
    return res
end

"""
Ordinary least squares for either FilteredUnknown[G|W]
"""
function fitcontiguouso(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * olspinv(extract(unk, chs), x, 1.0, lbls)
end
