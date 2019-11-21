using SparseArrays
using Polynomials
using LinearAlgebra
using TimerOutputs

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

    function TopHatFilter(ty::Type{<:FittingFilterType}, det::Detector, filt::AbstractMatrix{Float64}, wgts::AbstractVector{Float64})
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
        return new(ty, det, offsets, filts, wgts)
    end
end

Base.show(io::IO, thf::TopHatFilter) =
    print(io, "$(thf.filttype)[$(thf.detector)]")

"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=2.0)::TopHatFilter

Build the default top-hat filter for the specified detector with the specified top and base parameters.
"""
buildfilter(det::Detector, a::Float64 = 1.0, b::Float64 = 2.0, full::Bool = false)::TopHatFilter =
    buildfilter(VariableWidthFilter, det, a, b, full)

"""
    buildfilter(::Type{<:FittingFilterType}, det::Detector, a::Float64=1.0, b::Float64=1.0)::TopHatFilter

Build a top-hat-style filter for the specified detector with the specified top and base parameters. The
VariableWidthFilter and ConstantWidthFilter types are currently supported.
"""
function buildfilter(
    ty::Type{<:FittingFilterType},
    det::Detector,
    a::Float64 = 1.0, # Top
    b::Float64 = 2.0  # Base
)::TopHatFilter
    resf(::Type{VariableWidthFilter}, center, det) = resolution(center, det) # FWHM at center
    resf(::Type{ConstantWidthFilter}, center, det) = resolution(energy(n"Mn K-L3"), det) # FWHM at Mn Ka
    intersect(a, b, c, d) = max(0.0, min(b, d) - max(a, c)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1) * (-0.5 / (mina - minb)) + intersect(mina, maxa, e0, e1) / (maxa - mina) +
        intersect(maxa, maxb, e0, e1) * (-0.5 / (maxb - maxa))
    M(minb, mina) = channel(mina, det) - channel(minb, det)   # base length in channels
    N(mina, maxa) = (channel(maxa, det) - channel(mina, det))/2 # Top length in channels
    cc = channelcount(det)
    filt, wgts = zeros(Float64, (cc, cc)), zeros(Float64, cc)
    for ch1 in eachindex(wgts)
        center = 0.5 * (energy(ch1, det) + energy(ch1 + 1, det)) # midpoint of channel
        res = resf(ty, center, det)
        ea = ( center - 0.5 * a * res, center + 0.5 * a * res )
        eb = ( ea[1] - 0.5 * b * res, ea[2] + 0.5 * b * res )
        chmin, chmax = channel(eb[1], det), channel(eb[2], det)
        if (chmin>=1) && (chmax<=cc)
            for i in 0:(chmax-chmin)÷2
                filt[ch1, chmax-i] =
                    (filt[ch1, i+chmin] = filtint(energy(i+chmin, det), energy(i+chmin + 1, det), eb[1], ea[1], ea[2], eb[2]))
            end
            filt[ch1, chmax]= (filt[ch1, chmin]-=sum(filt[ch1,chmin:min(cc,chmax)])/2.0)
            @assert abs(sum(filt[ch1,:]))<1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i->filt[ch1,i]==filt[ch1,chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - O"
        end
        m, n = M(eb[1], ea[1]), N(ea...)
        @assert m>0.0 "m=$m!!!!"
        @assert n>0.0 "n=$n!!!!"
        wgts[ch1] = ((2.0*m+1.0)*(2.0*n))/(2.0*m+1.0+2.0*n) # Schamber's formula (see Statham 1977, Anal Chem 49, 14 )
    end
    @assert all(r->abs(sum(r))<1.0e-12, eachrow(filt)) "A filter does not sum to zero."
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
    b::Float64 = 5.0 # Width
)::TopHatFilter
    filtint(center, ee, gw) = exp(-0.5*((ee-center)/gw)^2)
    cc = channelcount(det)
    filt, wgts = zeros(Float64, (cc, cc)), zeros(Float64, cc)
    for ch1 in eachindex(wgts)
        center = energy(ch1, det) # midpoint of channel
        res = a * gaussianwidth(resolution(center, det))
        ex = ( center - 0.5 * b * res, center + 0.5 * b * res )
        chmin, chmax = channel(ex[1], det), channel(ex[2], det)
        if (chmin>=1) && (chmax<=cc)
            for i in 0:(chmax-chmin)÷2 # Ensure that it is symmetric
                filt[ch1, chmax-i] = (filt[ch1, chmin+i] = filtint(center, energy(chmin+i, det), res))
            end
            # Offset the Gaussian to ensure the sum is zero.
            filt[ch1, chmin:chmax] .-= sum(filt[ch1,chmin:chmax])/length(chmin:chmax)
            @assert abs(sum(filt[ch1,:]))<1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i->filt[ch1,i]==filt[ch1,chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - G"
        end
        wgts[ch1] = 1.0
    end
    return TopHatFilter(ty, det, filt, wgts)
end

"""
    FilteredLabel

An abstract type associated with labels of filtered spectrum data objects.
"""
abstract type FilteredLabel <: Label end

"""
    ReferenceLabel

A label associated with reference spectra.  The label encapsulates the original spectrum and the range of channels
represented by this reference object.
"""
abstract type ReferenceLabel <: FilteredLabel end

struct SpectrumLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
end
"""
    spectrum(fl::FilteredLabel)::Spectrum

The spectrum associated with a FilteredLabel-based type.
"""
spectrum(fl::FilteredLabel) = fl.spec

"""
   channels(rl::ReferenceLabel)::UnitRange{Int}

The range of channels associated with the specified ReferenceLabel.
"""
channels(rl::ReferenceLabel) = rl.roi

Base.show(io::IO, refLab::SpectrumLabel) = print(io::IO, "$(refLab.spec[:Name])[$(refLab.roi)]")
Base.isequal(rl1::SpectrumLabel, rl2::SpectrumLabel) = isequal(rl1.roi, rl2.roi) && isequal(rl1.spec, rl2.spec)

"""
    CharXRayLabel

A ReferenceLabel<:FilteredLabel  that Represents a reference spectrum associated with a set of characteristic x-rays
(CharXRay) objects over a contiguous range of spectrum channels.
"""
struct CharXRayLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
    xrays::Vector{CharXRay}
end

"""
   xrays(cl::CharXRayLabel)

A list of the X-rays associated with this CharXRayLabel.
"""
xrays(cl::CharXRayLabel) = cl.xrays

Base.show(io::IO, refLab::CharXRayLabel) = print(io::IO, "$(name(refLab.xrays))")
Base.isequal(rl1::CharXRayLabel, rl2::CharXRayLabel) =
    isequal(rl1.roi, rl2.roi) && isequal(rl1.xrays, rl2.xrays) && isequal(rl1.spec, rl2.spec)

"""
    EscapeLabel

A ReferenceLabel<:FilteredLabel that Represents a reference spectrum associated with an escape peak from a set of
characteristic x-rays (CharXRay) objects over a contiguous range of spectrum channels.
"""
struct EscapeLabel <: ReferenceLabel
    spec::Spectrum
    roi::UnitRange{Int}
    xrays::Vector{EscapeArtifact}
    EscapeLabel(spc::Spectrum, roi::UnitRange{Int}, escs::AbstractVector{EscapeArtifact}) =
        new(spc, roi, convert(Vector{EscapeArtifact},escs))
end

Base.show(io::IO, escl::EscapeLabel) = print(io, name(escl))
Base.isequal(el1::EscapeLabel, el2::EscapeLabel) =
    isequal(el1.roi, el2.roi) && isequal(el1.xrays, el2.xrays) && isequal(el1.spec, el2.spec)
NeXLCore.name(escl::EscapeLabel) = "Ecs[$(name([esc.xray for esc in escl.xrays]))]"



"""
    UnknownLabel

A FilteredLabel that represents the unknown spectrum.
"""
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
function _filter(
    spec::Spectrum,
    roi::UnitRange{Int},
    filter::TopHatFilter,
    tol::Float64
)
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
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
    return ( roiff, data[roiff], filtered[roiff] )
end

Base.convert(::AbstractVector{SpectrumFeature}, x::AbstractVector{CharXRay})::AbstractVector{SpectrumFeature} = x

charFeature(elm::Element, tr::Tuple{Vararg{Transition}}; minweight=1.0e-3, maxE=1.0e6) =
    convert(AbstractVector{SpectrumFeature},characteristic(elm, tr, minweight, maxE))

escapeFeature(elm::Element, trs::Tuple{Vararg{Transition}}; minweight=0.1, maxE=1.0e6, escape=n"Si K-L3") =
    convert(AbstractVector{SpectrumFeature}, map(tr->EscapeArtifact(tr,escape), characteristic(elm, trs, minweight, maxE)))

comptonFeature(elm::Element, trs::Tuple{Vararg{Transition}}; minweight=1.0e-3) =
    convert(AbstractVector{SpectrumFeature}, map(tr->ComptonArtifact(tr,escape), characteristic(elm, trs, minweight, maxE)))

function Base.filter(
    spec::Spectrum,
    det::Detector,
    cxrs::AbstractVector{SpectrumFeature},
    filter::TopHatFilter,
    scale::Float64=1.0;
    ampl::Float64=1.0e-4
)::Vector{FilteredReference}
    cxronly(lbl)::Vector{CharXRay} = Base.filter(l->l isa CharXRay, lbl)
    esconly(lbl)::Vector{EscapeArtifact} = Base.filter(l->l isa EscapeArtifact, lbl)
    buildlabel(spec, roi, lbl)::ReferenceLabel =
        if all(l->l isa EscapeArtifact, lbl)
            return EscapeLabel(spec, roi, esconly(lbl))
        else
            return CharXRayLabel(spec, roi, cxronly(lbl))
        end
    mb(spec, roi, lbl) =
        all(l->l isa EscapeArtifact, lbl) ? #
            modelBackground(spec, roi) : #
            modelBackground(spec, roi, inner(brightest(cxronly(lbl))))
    res=[]
    for (lbl, roi) in labeledextents(cxrs, det, ampl)
        charonly = spec[roi] - mb(spec, roi, lbl)
        f = _filter(spec, roi, filter, 1.0e-6)
        push!(res, FilteredReference(buildlabel(spec, roi, lbl), scale, roi, f..., charonly, filter.weights[(roi.start+roi.stop)÷2]))
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
    charonly = spec[roi] - modelBackground(spec, roi, ashell)
    f = _filter(spec, roi, filter, tol)
    return FilteredReference(SpectrumLabel(spec, roi), scale, roi, f..., charonly, filter.weights[(roi.start+roi.stop)÷2])
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
    charonly = spec[roi] - modelBackground(spec, roi)
    f = _filter(spec, roi, filter, tol)
    return FilteredReference(SpectrumLabel(spec, roi), scale, roi, f..., charonly, filter.weights[(roi.start+roi.stop)÷2])
end

"""
    FilteredUnknown

A FilteredDatum representing the unknown.  There are two types of FilteredUnknown - FilteredUnknownG and
FilteredUnknownW for "generalized" and "weighted" model unknowns.  The distinction is desirable since the
full covariance calculation for the generalized model is so expensive relative to the weighted approximation.
"""
abstract type FilteredUnknown <: FilteredDatum end


"""
    FilteredUnknownG

Represents the unknown in a filter fit using the full generalized fitting model.  This model is expensive to
calculate but uses the full generalized linear fitting model which produces the correct fit covariances.
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
    function covariance(roi, data, filter) # Calculated the full covariance matrix
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
    FilteredUnknownW

Represents the unknown in a filter fit using the weighted fitting model.  This is an approximation that produces over
optimistic resulting covariance matrix.
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
    range(i) = filter.offsets[i] : filter.offsets[i] + length(filter.filters[i]) - 1
    apply(filter::TopHatFilter, data::AbstractVector{Float64}) =
        [ dot(filter.filters[i], view(data, range(i))) for i in eachindex(data) ]
    covariance(filter, data) = # The diagnonal elements of the full covariance matrix
        [ dot(filter.filters[r] .* view(data,range(r)), filter.filters[r]) for r in eachindex(data) ]
    data = counts(spec, Float64, true)
    # Compute the filtered data
    filtered = apply(filter, data)
    roi = eachindex(filtered)
    # max(d,1.0) is necessary to ensure the variances are positive
    covar = covariance(filter, map(d -> max(d, 1.0), data))
    return FilteredUnknownW(UnknownLabel(spec), scale, roi, roi, data, filtered, covar)
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
    covariance(fd::FilteredUnknownG, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.  <code>roi</code> must
fully encompass the filtered edata in <code>fd</code>.
"""
function covariance(fd::FilteredUnknownG, roi::UnitRange{Int})
    # Unknown case
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return Matrix(fd.covariance[nz,nz])
end

"""
    covariance(fd::FilteredUnknownW, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance diagnonal elements over the specified range of channels.
<code>roi</code> must fully encompass the filtered edata in <code>fd</code>.
"""
function covariance(fd::FilteredUnknownW, roi::UnitRange{Int})
    # Unknown case
    nz = roi.start-fd.ffroi.start+1:roi.stop-fd.ffroi.start+1
    return fd.covariance[nz]
end

# _buildXXX - Helpers for fitcontiguousX functions...
function _buildmodel(ffs::Array{FilteredReference}, chs::UnitRange{Int})::Matrix{Float64}
    x = Matrix{Float64}(undef, (length(chs), length(ffs)))
    for i in eachindex(ffs)
        x[:, i] = extract(ffs[i], chs)
    end
    return x
end
_buildlabels(ffs::Array{FilteredReference}) = collect(ff.identifier for ff in ffs)
_buildscale(unk::FilteredUnknown, ffs::Array{FilteredReference}) =
     Diagonal( [ unk.scale / ffs[i].scale for i in eachindex(ffs) ] )

"""
Generalized least squares (my implementation)
"""
function fitcontiguousg(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glssvd(extract(unk, chs), x, covariance(unk, chs), lbls)
end

"""
Generalized least squares (pseudo-inverse)
"""
function fitcontiguousp(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, covariance(unk, chs), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * glspinv(extract(unk, chs), x, Diagonal(diag(covariance(unk, chs))), lbls)
end

"""
Weighted least squares using the diagonal from the FilteredUnknownG covariance.
"""
function fitcontiguousw2(unk::FilteredUnknownG, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * wlssvd(extract(unk, chs), x, diag(covariance(unk, chs)), lbls)
end

"""
Weighted least squares for FilteredUnknownW
"""
function fitcontiguousww(unk::FilteredUnknownW, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
    wgts = [ ff.covscale for ff in ffs ]
    return scale * wlspinv(extract(unk, chs), x, covariance(unk, chs), wgts, lbls)
end

"""
Ordinary least squares for either FilteredUnknown[G|W]
"""
function fitcontiguouso(unk::FilteredUnknown, ffs::Array{FilteredReference}, chs::UnitRange{Int})::UncertainValues
    x, lbls, scale = _buildmodel(ffs,chs), _buildlabels(ffs), _buildscale(unk, ffs)
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
    peakback::Dict{<:ReferenceLabel,NTuple{2, Float64}}
end

"""
    kratios(ffr::FilterFitResult)

The k-ratios as a UncertainValues object
"""
kratios(ffr::FilterFitResult)::UncertainValues = ffr.kratios
unknown(ffr::FilterFitResult)::UnknownLabel = ffr.label
residual(ffr::FilterFitResult)::Vector{Float64} = ffr.residual

NeXLCore.asa(::Type{DataFrame}, ffrs::Array{FilterFitResult}, withUnc=false)::DataFrame =
    asa(DataFrame, [r.kratios for r in ffrs], withUnc)

function NeXLCore.asa(::Type{DataFrame}, ffr::FilterFitResult)::DataFrame
    lbl, klbl, kr, dkr, roi1, roi2, peak, back =
        UnknownLabel[], ReferenceLabel[], Float64[], Float64[], Int[], Int[], Float64[], Float64[]
    for kl in labels(ffr.kratios)
        push!(lbl,ffr.label)
        push!(klbl, kl)
        push!(roi1, kl.roi.start)
        push!(roi2, kl.roi.stop)
        push!(kr, value(kl, ffr.kratios))
        push!(dkr, σ(kl, ffr.kratios))
        pb = ffr.peakback[kl]
        push!(peak, pb[1])
        push!(back, pb[2])
    end
    return DataFrame(Label=lbl, Feature=klbl, Start=roi1, Stop=roi2, K=kr, dK=dkr, Peak=peak, Back=back)
end

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
)::Dict{<:ReferenceLabel,NTuple{2, Float64}}
    res = Dict{ReferenceLabel,NTuple{2, Float64}}()
    for ff in ffs
        su = sum(unk.data[ff.roi])
        res[ff.identifier] = (su, su - sum((value(ff.identifier, kr) * ff.scale / unk.scale) * ff.charonly))
    end
    return res
end

"""
    filterfit(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg = fitcontiguousp)::FilterFitResult
    fitrois = ascontiguous(map(fd->fd.ffroi, ffs))  # Divide up into contiguous rois
    @info "Fitting $(length(ffs)) references in $(length(fitrois)) blocks - $fitrois"
    kr = cat(map(fr->alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, ffs), fr), fitrois)) # fit each roi
    resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, resid, pb) # Return results
end

"""
    filterfit_all(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues


Filter fit the unknown against ffs, an array of FilteredReference and return the result as an FilterFitResult object.
Fits all the references as one large linear equation rather than breaking them up into contiguous blocks.
By default use the generalized LLSQ fitting (pseudo-inverse implementation).
"""
function filterfit_all(unk::FilteredUnknownG, ffs::Array{FilteredReference}, alg = fitcontiguousp)::FilterFitResult
    fr = minimum(ff.ffroi.start for ff in ffs)-20:maximum(ff.ffroi.stop for ff in ffs)+20  # Make one large roi
    @info "Fitting $(length(ffs)) references in 1 block - $fr"
    kr = alg(unk, ffs, fr) # Fit the roi
    resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, resid, pb) # Return results
end


"""
    filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg=fitcontiguousw)::UncertainValues

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
function filterfit(unk::FilteredUnknownW, ffs::Array{FilteredReference}, alg = fitcontiguousww, forcezeros=true)::FilterFitResult
    trimmed, refit, removed, retained = copy(ffs), true, Vector{UncertainValues}(), nothing # start with all the FilteredReference
    while refit
        fitrois = ascontiguous(map(fd->fd.ffroi, trimmed))
        # @info "Fitting $(length(trimmed)) references in $(length(fitrois)) blocks - $fitrois"
        retained = map(fr->alg(unk, filter(ff -> length(intersect(fr, ff.ffroi)) > 0, trimmed), fr), fitrois)
        kr = cat(retained)
        refit = false
        for lbl in labels(kr)
            if value(lbl, kr) <= 0.0
                splice!(trimmed, findfirst(ff->ff.identifier==lbl, trimmed))
                push!(removed, uvs([lbl],[forcezeros ? 0.0 : value(lbl, kr)],reshape([σ(lbl, kr)],(1,1))))
                refit=true
            end
        end
    end # while
    kr = cat(append!(retained, removed))
    resid, pb = _computeResidual(unk, ffs, kr), _computecounts(unk, ffs, kr)
    return FilterFitResult(unk.identifier, kr, unk.roi, unk.data, resid, pb)
end

"""
    filteredresidual(fit::UncertainValues, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
function filteredresidual(fit::FilterFitResult, unk::FilteredUnknown, ffs::Array{FilteredReference})::Vector{Float64}
    scaled(ff) = (value(ff.identifier, fit) * (ff.scale / unk.scale)) * extract(ff, unk.ffroi)
    return unk.filtered - mapreduce(scaled, +, ffs)
end

fit(ty::Type{FilteredUnknownW}, unk::Spectrum, filt::TopHatFilter, refs::Vector{FilteredReference}, forcezeros=true) =
    filterfit(filter(ty, unk, filt, 1.0/dose(unk)), refs, fitcontiguousww, forcezeros)

fit(ty::Type{FilteredUnknownG}, unk::Spectrum, filt::TopHatFilter, refs::Vector{FilteredReference}, forcezeros=true) =
    filterfit(filter(ty, unk, filt, 1.0/dose(unk)), refs, fitcontiguousw, forcezeros)

fit(unk::Spectrum, filt::TopHatFilter, refs::Vector{FilteredReference}, forcezeros=true) =
    fit(FilteredUnknownW, unk, filt, refs, forcezeros)
