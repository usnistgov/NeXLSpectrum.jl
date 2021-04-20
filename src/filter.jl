
"""
    TopHatFilterType

Represents different fitting filter models.  A fitting filter is a linear operator that is 1) symmetric about the
center, and 2) the sum of the elements must be zero.  The standard filter is a top-hat shape, although Gaussian,
triangular, Savitsky-Golay and other shapes are also possible.
"""
abstract type TopHatFilterType end

"""
    VariableWidthFilter

A top-hat filter that varies in width with the FWHM of the detector.
"""
struct VariableWidthFilter <: TopHatFilterType end # Variable width filter

"""
    ConstantWidthFilter

A top-hat filter that has constant width determined by FWHM at Mn Kα for all channels.
"""
struct ConstantWidthFilter <: TopHatFilterType end # Constant width filter

"""
    GaussianFilter

A Gaussian-shaped filter that varies in width with the FWHM of the detector.  The Gaussian is offset to ensure the
sum of the filter elements is zero.
"""
struct GaussianFilter <: TopHatFilterType end # A variable width Gaussian filter

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

Notes on memory and code optimization:
The filter matrix is banded diagonal.  Approximately, 2.5% of the elements are non-zero.  This suggest use of the
BandedMatrix type.  The most expensive operation is calculating F⋅D⋅Fᵀ, the covariance matrix of the filtered data.
D is a diagonal matrix and so computing each element in F⋅D⋅Fᵀ reduces to a sum over a single variable.
Furthermore, the weighted least squares fit doesn't require the full F⋅D⋅Fᵀ, just diag(F⋅D⋅Fᵀ).  However, it turns out
that we can do better implementing our own banded matrix type largely because D is fully diagonal and the matrix
product F⋅D⋅Fᵀ reduces down to a sum over a single variable.  The product F⋅d and F⋅D⋅Fᵀ are readily
implemented as element-by-element multiplies and sums.  Thus storing the filter as offsets and row filters is
efficient in both memory and CPU use.
"""
struct TopHatFilter
    filttype::Type{<:TopHatFilterType}
    detector::Detector
    offsets::Vector{Int}  # offset to start of filter in row
    filters::Vector{Vector{Float64}} # Filter data
    weights::AbstractVector{Float64} # Correlation compensation weights

    function TopHatFilter(
        ty::Type{<:TopHatFilterType},
        det::Detector,
        filt::AbstractMatrix{Float64},
        wgts::AbstractVector{Float64},
    )
        # Extract the contiguous non-zero elements
        dim = size(filt)[1]
        offsets = zeros(Int, dim)
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

function filterdata(filt::TopHatFilter, region::AbstractUnitRange{Int})::Matrix{Float64}
    res = zeros(Float64, length(region), length(region))
    foreach(r -> res[r-first(region)+1, :] = filterdata(filt, r)[region], region)
    return res
end

filterdata(filt::TopHatFilter) = #
    filterdata(filt, eachindex(filt.filters))

"""
    filteredcovar(filt::TopHatFilter, specdata::Vector{Float64}, row::Int, col::Int)::Float64

Compute the covariance matrix entry for the specified row and column.  The resulting covariance matrix is equal
to F⋅Ω⋅Fᵀ where F is the filter and Ω is the covariance matrix of the spectrum data.  Since each channel in the spectrum
is 1) independent; 2) Poisson distributed, Ωᵢⱼ = specdata[i] if i==j, 0 otherwise. Because Ω is diagonal, the matrix
multiplication reduces to sum over a single index.  Often this sum is zero because the non-zero portions of the
`row` and `col` filters don't intersect.

Note:  `specdata` should be preprocessed so that no element is less than or equal to zero.
"""
function filteredcovar(
    filt::TopHatFilter,
    specdata::Vector{Float64},
    i::Int,
    l::Int,
)::Float64
    function dot3(a, b, c) # vmapreduce((ai,bi,ci)->ai*bi*ci, +, a, b, c)
        sum = zero(eltype(a))
        @avx for i in eachindex(a)  # @avx takes overall time from 670 μs down to 430 μs
            sum += a[i] * b[i] * c[i]
        end
        return sum
    end
    fi, fl, oi, ol = filt.filters[i], filt.filters[l], filt.offsets[i], filt.offsets[l]
    roi = max(oi, ol):min(oi + length(fi) - 1, ol + length(fl) - 1) # The ROI over which both filters are non-zero.
    return length(roi) > 0 ? #
           dot3(
        view(fi, roi.start-oi+1:roi.stop-oi+1),
        view(specdata, roi),
        view(fl, roi.start-ol+1:roi.stop-ol+1),
    ) : 0.0
end

"""
    filtereddatum(filt::TopHatFilter, specdata::Vector{Float64}, ch::Int)::Float64 =

Compute a single channel in the filtered spectrum.
"""
filtereddatum(filt::TopHatFilter, specdata::Vector{Float64}, i::Int)::Float64 = dot(
    filt.filters[i],
    view(specdata, filt.offsets[i]:filt.offsets[i]+length(filt.filters[i])-1),
)


Base.size(filt::TopHatFilter) = (length(filt.filters), length(filt.filters))
Base.length(filt::TopHatFilter) = length(filt.filters)

Base.show(io::IO, thf::TopHatFilter) = print(io, "$(thf.filttype)[$(thf.detector)]")

"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=2.0)::TopHatFilter

Build the default top-hat filter for the specified detector with the specified top and base parameters.
"""
buildfilter(det::Detector, a::Float64 = 1.0, b::Float64 = 2.0)::TopHatFilter =
    buildfilter(VariableWidthFilter, det, a, b)

"""
    buildfilter(::Type{<:TopHatFilterType}, det::Detector, a::Float64=1.0, b::Float64=1.0)::TopHatFilter

Build a top-hat-style filter for the specified detector with the specified top and base parameters. The
VariableWidthFilter and ConstantWidthFilter types are currently supported.
"""
function buildfilter(
    ty::Type{<:TopHatFilterType},
    det::Detector,
    a::Float64 = 1.0, # Top
    b::Float64 = 2.0,  # Base
)::TopHatFilter
    resf(::Type{VariableWidthFilter}, center, det) = resolution(center, det) # FWHM at center
    resf(::Type{ConstantWidthFilter}, center, det) = resolution(energy(n"Mn K-L3"), det) # FWHM at Mn Ka
    intersect(aa, bb, cc, dd) = max(0.0, min(bb, dd) - max(aa, cc)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1) * (-0.5 / (mina - minb)) +
        intersect(mina, maxa, e0, e1) / (maxa - mina) +
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
                filt[ch1, chmax-i] = (
                    filt[ch1, i+chmin] = filtint(
                        energy(i + chmin, det),
                        energy(i + chmin + 1, det),
                        eb[1],
                        ea[1],
                        ea[2],
                        eb[2],
                    )
                )
            end
            filt[ch1, chmax] =
                (filt[ch1, chmin] -= sum(filt[ch1, chmin:min(cc, chmax)]) / 2.0)
            @assert abs(sum(filt[ch1, :])) < 1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i -> filt[ch1, i] == filt[ch1, chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - O"
        end
        m, n = M(ea[2], ea[1]), N(eb[1], ea[1])
        @assert m > 0.0 "m=$(m)!!!!"
        @assert n > 0.0 "n=$(n)!!!!"
        wgts[ch1] = 2.87 + 1.758e-4 * energy(ch1, det)
        # wgts[ch1] = 2.0*n*m/(n+2.0*m)  # Schamber's formula (see Schamber 1997, note the formula in Statham 1977, Anal Chem 49, 14 doesn't seem to work.)
    end
    @assert all(r -> abs(sum(r)) < 1.0e-12, eachrow(filt)) "A filter does not sum to zero."
    return TopHatFilter(ty, det, filt, wgts)
end

"""
    buildfilter(::Type{GaussianFilter}, det::Detector, a::Float64=1.0, b::Float64=5.0)::TopHatFilter

Build a top-hat filter with Gaussian shape whose width varies with the detector's resolution as a function of X-ray
energy for the specified detector with the specified top and base parameters. The `a` parameter corresponds
to the filter width relative to the detector resolution expressed as Gaussian width.  So `a=1` is a filter
whose width equals the detector resolution at each energy.  The `b` parameter is the extent of the
filter in Gaussian widths.  The default `a=1, b=5` corresponds to a  filter that has the same resolution
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
                filt[ch1, chmax-i] =
                    (filt[ch1, chmin+i] = filtint(center, energy(chmin + i, det), res))
            end
            # Offset the Gaussian to ensure the sum is zero.
            filt[ch1, chmin:chmax] .-= sum(filt[ch1, chmin:chmax]) / length(chmin:chmax)
            @assert abs(sum(filt[ch1, :])) < 1.0e-12 "Filter $ch1 does not sum to zero."
            @assert all(i -> filt[ch1, i] == filt[ch1, chmax-(i-chmin)], chmin:chmax) "The $ch1-th filter is not symmetric - G"
        end
        wgts[ch1] = 2.87 + 1.758e-4 * energy(ch1, det)
    end
    @info "The uncertainty estimates will be about a factor of three low for the Gaussian filter."
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
    label::ReferenceLabel # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    charonly::Vector{Float64} # Background corrected intensity data over roi
    covscale::Float64  # Schamber's variance scale factor
end

Base.show(io::IO, fd::FilteredReference) = print(io, "Reference[$(fd.label)]")

"""
    _filter(
      spec::Spectrum,
      roi::UnitRange{Int},
      thf::TopHatFilter,
      tol::Float64
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.
"""
function _filter(spec::Spectrum, roi::UnitRange{Int}, thf::TopHatFilter, tol::Float64)
    # Extract the spectrum data as Float64 to match the filter
    data = counts(spec, Float64, true)
    # Determine tangents to the two background end points
    tangents = map(st -> estimatebackground(data, st, 5, 2), (roi.start, roi.stop))
    # Replace the non-ROI channels with extensions of the tangent functions
    data[1:roi.start-1] = map(tangents[1], 1-roi.start:-1)
    data[roi.stop+1:end] = map(tangents[2], 1:length(data)-roi.stop)
    # Compute the filtered data
    filtered = [filtereddatum(thf, data, i) for i in eachindex(data)]
    maxval = maximum(filtered)
    ff = findfirst(f -> abs(f) > tol * maxval, filtered)
    if isnothing(ff)
        error("There doesn't appear to be any data in $(spec[:Name]))")
    end
    roiff = ff:findlast(f -> abs(f) > tol * maxval, filtered)
    return (roiff, data[roiff], filtered[roiff])
end

function tophatfilter(
    escLabel::EscapeLabel,
    thf::TopHatFilter,
    scale::Float64,
    tol::Float64 = 1.0e-6,
)::FilteredReference
    spec, roi = escLabel.spectrum, escLabel.roi
    charonly = spec[roi] - modelBackground(spec, roi)
    f = _filter(spec, roi, thf, tol)
    return FilteredReference(
        escLabel,
        scale,
        roi,
        f...,
        charonly,
        thf.weights[(roi.start+roi.stop)÷2],
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
labeledextents(elm::Element, det::Detector, ampl::Float64, maxE::Float64 = 1.0e6) =
    labeledextents(
        isvisible(characteristic(elm, alltransitions, ampl, maxE), det),
        det,
        ampl,
    )

"""
    charXRayLabels(#
      spec::Spectrum, #
      elm::Element, #
      allElms::AbstractVector{Element}, #
      det::Detector, #
      ampl::Float64, #
      maxE::Float64=1.0e6)::Vector{SpectrumFeature}

Creates a vector of CharXRayLabel objects associated with 'elm' for a spectrum containing the elements
'allElms' assuming that it was collected on 'det'.  ROIs in which other elements from 'allElms'
interfere with 'elm' will not be included.
"""
function charXRayLabels(#
    spec::Spectrum, #
    elm::Element, #
    allElms::AbstractVector{Element}, #
    det::Detector, #
    maxE::Float64 = 1.0e6;
    ampl::Float64 = 1.0e-5, #
)::Vector{ReferenceLabel}
    suitable = suitablefor(elm, allElms, det, maxE, ampl)
    length(suitable) > 0 || @error "The spectrum $(name(spec)) provides no suitable ROIs for the element $(symbol(elm))."
    # Find all the ROIs associated with other elements
    return [ CharXRayLabel(spec, roi, xrays) for (xrays, roi) in suitable ]
end

"""
    suitablefor(elm::Element, allElms::AbstractVector{Element}, det::Detector, maxE::Float64 = 1.0e6, ampl::Float64 = 1.0e-5)::Vector{Tuple{Vector{CharXRay}, UnitRange{Int64}}}

Given a material containing `allElms` which ROIs for element `elm` is the material a suitable reference on the detector `det`?  
This function provides informational, warning and error messages depending upon the suitability of the material. 
"""
function suitablefor( #
    elm::Element, #
    allElms::AbstractVector{Element}, #
    det::Detector, #
    maxE::Float64 = 1.0e6,
    ampl::Float64 = 1.0e-5, #
    warnme::Bool = true #
)::Vector{Tuple{Vector{CharXRay}, UnitRange{Int64}}}
    # Find all the ROIs associated with other elements
    if elm in allElms 
        otherElms = filter(a -> a ≠ elm, allElms)
        res = if length(otherElms) > 0
            urs = mapreduce(ae -> extents(ae, det, ampl), append!, otherElms)
            # Find elm's ROIs that don't intersect another element's ROI
            filter(labeledextents(elm, det, ampl, maxE)) do lx 
                furs = filter(ur -> length(intersect(ur, lx[2])) > 0, urs)
                warnme && (!isempty(furs)) && @info "A material containing $(join(symbol.(allElms),", "," and ")) is not a suitable reference for \"$(lx[1])\" due to $(length(furs)==1 ? "a peak interference." : "$(length(furs)) peak interferences.")"
                isempty(furs)
            end
        else
            labeledextents(elm, det, ampl, maxE)
        end
        length(res) == 0 && warnme && @warn "A material containing $(join(symbol.(allElms),", "," and ")) provides no references for $(symbol(elm))."
    else
        @error "The element $(symbol(elm)) is not contained within the elements $(join(symbol.(allElms),", "," and "))."
        res = Tuple{Vector{CharXRay}, UnitRange{Int64}}[]
    end
    return res
end
function suitablefor( #
    elm::Element, #
    mat::Material, #
    det::Detector, #
    maxE::Float64 = 1.0e6,
    ampl::Float64 = 1.0e-5, #
    warnme::Bool = true #
)
    suitablefor(elm, collect(keys(mat)), det, maxE, ampl, warnme)
end

"""
    suitability(elm::Element, mats::AbstractSet{<:Material}, det::Detector; maxE=30.0e3)
    suitability(elm::Element, det::Detector; maxE=30.0e3, minC=0.1)

Tabulates the characteristic X-ray peaks for the `Element` for which there are suitable materials 
in `mats` for the specified detector.  The second form uses a default set of Materials in the
file NeXLCore "standards.txt".

This function is helpful for determining which `Material`s are suitable to act as 
fitting standards for the specified Element.  It shows how NeXLSpectrum will break up the 
characteristic peaks associated with `elm` into contiguous regions each of which will be 
fit independently. NeXLSpectrum attempts to break each element into as many independent 
regions as possible dependent on the resolution of the specified `EDSDetector`.  If there
is an interference between one of the other elements in the `Material` and `elm` then
this peak will not be suitable as a fitting standard.  However, it can be used as a 
similar standard.
"""
function suitability(elm::Element, mats::AbstractSet{<:Material}, det::EDSDetector; maxE=30.0e3, latex=false)
    wrap(matname, ltx) = ltx ? "\\ce{$matname}" : matname
    ex, chk = latex ? ( "\\xmark", "\\checkmark" ) : ("✗", "✓")
    mats = filter(m->value(m[elm])>0.0, mats)
    # Find all elemental ROIs
    rois = Dict{Vector{CharXRay}, Vector{Material}}()
    for (cxrs, ur) in NeXLSpectrum.suitablefor(elm, pure(elm), det, maxE, 1.0e-5, false)
      rois[cxrs] = Material[]
    end
    for m in mats
      for (cxrs, ur) in NeXLSpectrum.suitablefor(elm, m, det, maxE, 1.0e-5, false)
        push!(rois[cxrs], m) 
      end
    end
    cxrss = collect(keys(rois))
    sort!(cxrss, lt = (a,b) -> energy(brightest(a))<energy(brightest(b)))
    res = DataFrame(:Material=>String[], :Count=>Int[], ( Symbol(cxrs)=>String[] for cxrs in cxrss )...)
    for m in mats
      push!(res, [wrap(name(m), latex), count(cxrs->m in rois[cxrs], cxrss), (m in rois[cxrs] ? chk : ex for cxrs in cxrss)...])
    end
    sort!(res, :Count, rev=true)
end
function suitability(elm::Element, det::EDSDetector; maxE=30.0e3, minC=0.1, latex=false)
    suitability(elm, getstandards(elm, minC), det, maxE=maxE, latex=latex)
end

function escapeextents(elm::Element, det::Detector, ampl::Float64, maxE::Float64)
    vis = isvisible(characteristic(elm, alltransitions, ampl, maxE), det)
    escapeextents(vis, det, ampl, maxE)
end


"""
    escapeLabels(#
      spec::Spectrum, #
      elm::Element, #
      allElms::AbstractVector{Element}, #
      det::Detector, #
      ampl::Float64, #
      maxE::Float64=1.0e6
    )::Vector{EscapeLabel}

Creates a vector EscapeLabel objects associated with 'elm' for a spectrum containing the elements
'allElms' assuming that it was collected on 'det'.  ROIs in which other elements from 'allElms'
interfere with 'elm' will not be included.
"""
function escapeLabels(#
    spec::Spectrum, #
    elm::Element, #
    allElms::AbstractVector{Element}, #
    det::Detector, #
    maxE::Float64 = 1.0e6;
    ampl::Float64 = 1.0e-5, #
)::Vector{ReferenceLabel}
    @assert elm in allElms "$elm must be in $allElms."
    intersects(urs, roi) = !isnothing(findfirst(ur -> length(intersect(ur, roi)) > 0, urs))
    # Find all the ROIs associated with the elements
    urs = collect(Iterators.flatten(extents(ae, det, ampl) for ae in allElms))
    # Find elm's ROIs that don't intersect an element's primary ROI
    lxs = filter(lx -> !intersects(urs, lx[2]), escapeextents(elm, det, ampl, maxE))
    isempty(lxs) || @warn "There are no escape peaks for $elm that don't interfere with one or more of $allElms."
    return ReferenceLabel[EscapeLabel(spec, roi, xrays) for (xrays, roi) in lxs]
end

"""
    tophatfilter(
        charLabel::CharXRayLabel,
        thf::TopHatFilter,
        scale::Float64 = 1.0,
        tol::Float64 = 1.0e-6,
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.  Use a simple
edge-based background model.
"""
function tophatfilter(
    charLabel::CharXRayLabel,
    filt::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::FilteredReference
    ashellof(xrays) = inner(xrays[argmax(jumpratio.(inner.(xrays)))])
    spec, roi, ashell = charLabel.spectrum, charLabel.roi, ashellof(charLabel.xrays)
    return FilteredReference(
        charLabel,
        scale,
        roi,
        _filter(spec, roi, filt, tol)...,
        spec[roi] - modelBackground(spec, roi, ashell),
        filt.weights[(roi.start+roi.stop)÷2],
    )
end

function tophatfilter(
    labels::AbstractVector{ReferenceLabel},
    filt::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::Vector{FilteredReference} 
    lbls = filter(labels) do lbl
        res = (lbl.roi.start >= 1) && (lbl.roi.stop <= length(filt)) 
        (!res) && @warn "The ROI $(lbl.roi) not fully contained on [1, $(length(filt))] for $lbl."
        w = (weight(brightest(lbl.xrays)) > 1.0e-3)
        (!w) && @warn "No sufficiently bright X-rays for $lbl."
        res && w
    end    
    return [ tophatfilter(lbl, filt, scale, tol) for lbl in lbls ]
end

"""
    tophatfilter(
      reflabel::ReferenceLabel,
      roi::UnitRange{Int},
      thf::TopHatFilter,
      scale = 1.0,
      tol = 1.0e-6
    )::FilteredReference

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter. Use a naive
linear background model.
"""
function tophatfilter(
    reflabel::ReferenceLabel,
    thf::TopHatFilter,
    scale::Float64 = 1.0,
    tol::Float64 = 1.0e-6,
)::FilteredReference
    spec, roi = reflabel.spectrum, reflabel.roi
    return FilteredReference(
        reflabel,
        scale,
        roi,
        _filter(spec, roi, thf, tol)...,
        spec[roi] - modelBackground(spec, roi),
        thf.weights[(roi.start+roi.stop)÷2],
    )
end

"""
    FilteredUnknown

A FilteredDatum representing the unknown.  There are two types of FilteredUnknown - FilteredUnknownG and
FilteredUnknownW for "generalized" and "weighted" model unknowns.  The distinction is desirable since the
full covariance calculation for the generalized model is so expensive relative to the weighted approximation.
"""
abstract type FilteredUnknown <: FilteredDatum end

Base.show(io::IO, fd::FilteredUnknown) = print(io, fd.label)

"""
    extract(fd::FilteredReference, roi::UnitRange{Int})

Extract the filtered data representing the specified range.  `roi` must fully encompass the filtered
data in `fd`.
"""
function NeXLUncertainties.extract(fd::FilteredReference, roi::UnitRange{Int})
    @assert fd.ffroi.start >= roi.start "$(fd.ffroi.start) < $(roi.start) in $(fd)"
    @assert fd.ffroi.stop <= roi.stop "$(fd.ffroi.stop) > $(roi.stop) in $(fd)"
    data = zeros(Float64, length(roi))
    nz = fd.ffroi.start-roi.start+1:fd.ffroi.stop-roi.start+1
    data[nz] = fd.filtered
    return data
end

"""
    extract(fd::FilteredUnknown, roi::UnitRange{Int})::AbstractVector{Float64}

Extract the filtered data representing the specified range.  `roi` must be fully contained within the
filtered data in `fd`.
"""
NeXLUncertainties.extract(
    fd::FilteredUnknown,
    roi::UnitRange{Int},
)::AbstractVector{Float64} = view(fd.filtered, roi)

_buildlabels(ffs::AbstractVector{FilteredReference}) = collect(ff.label for ff in ffs)
_buildscale(unk::FilteredUnknown, ffs::AbstractVector{FilteredReference}) =
    Diagonal([unk.scale / ffs[i].scale for i in eachindex(ffs)])

# Internal: Computes the residual spectrum based on the fit k-ratios
function _computeResidual(
    unk::FilteredUnknown,
    ffs::AbstractVector{FilteredReference},
    kr::UncertainValues,
)
    res = copy(unk.data)
    for ff in ffs
        res[ff.roi] -=
            (NeXLUncertainties.value(ff.label, kr) * ff.scale / unk.scale) * ff.charonly
    end
    return res
end

# Internal: Computes the peak and background count based on the fit k-ratios
function _computecounts( #
    unk::FilteredUnknown, #
    ffs::AbstractVector{FilteredReference}, #
    kr::UncertainValues, #
)::Dict{<:ReferenceLabel,NTuple{2,Float64}}
    res = Dict{ReferenceLabel,NTuple{2,Float64}}()
    for ff in ffs
        su = sum(unk.data[ff.roi])
        res[ff.label] = (
            su,
            su - sum(
                (NeXLUncertainties.value(ff.label, kr) * ff.scale / unk.scale) *
                ff.charonly,
            ),
        )
    end
    return res
end

"""
Ordinary least squares for either FilteredUnknown[G|W]
"""
function fitcontiguouso(
    unk::FilteredUnknown,
    ffs::AbstractVector{FilteredReference},
    chs::UnitRange{Int},
)::UncertainValues
    x, lbls, scale = _buildmodel(ffs, chs), _buildlabels(ffs), _buildscale(unk, ffs)
    return scale * olspinv(extract(unk, chs), x, 1.0, lbls)
end

"""
    selectBestReferences(refs::AbstractVector{FilteredReference})::Vector{FilteredReference}

For each roi in each element, pick the best FilteredReference with the highest
intensity in the characteristic-only data.
"""
function selectBestReferences(
    refs::AbstractVector{FilteredReference},
)::Vector{FilteredReference}
    rois = Dict{UnitRange,Tuple{FilteredReference,Float64}}()
    elms = unique(element(ref.label) for ref in refs)
    for elm in elms
        for ref in filter(ref -> element(ref.label) == elm, refs)
            comp = get(ref.label.spectrum, :Composition, missing)
            # Pick the reference with the largest charonly value but prefer pure elements over compounds
            mrf = maximum(ref.charonly) * (ismissing(comp) ? 0.01 : normalized(comp, elm)^2)
            if (!haskey(rois, ref.roi)) || (mrf > rois[ref.roi][2])
                rois[ref.roi] = (ref, mrf)
            end
        end
    end
    return [v[1] for v in values(rois)]
end

function filterreference(
    filt::TopHatFilter,
    spec::Spectrum,
    elm::Element,
    elms::AbstractArray{Element};
    props::Dict{Symbol,<:Any} = Dict{Symbol,Any}(),
    withEsc::Bool = false,
)::Vector{FilteredReference}
    @assert elm in elms
    cprops = merge(spec.properties, props)
    if !haskey(cprops, :Detector)
        cprops[:Detector] = filt.detector
    end
    # Creates a list of unobstructed ROIs for elm as CharXRayLabel objects
    lbls = charXRayLabels(spec, elm, elms, cprops[:Detector], cprops[:BeamEnergy])
    if withEsc
        append!(lbls, escapeLabels(spec, elm, elms, cprops[:Detector], cprops[:BeamEnergy]))
    end
    # Filters the spectrum and returns a list of FilteredReference objects, one per ROI
    return tophatfilter(lbls, filt, 1.0 / (cprops[:ProbeCurrent] * cprops[:LiveTime]))
end

function filterreference(
    filt::TopHatFilter,
    spec::Spectrum,
    elm::Element,
    comp::Material;
    props::Dict{Symbol,<:Any} = Dict{Symbol,Any}(),
    withEsc::Bool = false,
)::Vector{FilteredReference}
    spec[:Composition] = get(spec, :Composition, comp)
    filterreference(filt, spec, elm, collect(keys(comp)), props = props, withEsc = withEsc)
end

filterreferences(
    filt::TopHatFilter,
    refs::Tuple{Spectrum,Element,Material}...;
    props::Dict{Symbol,<:Any} = Dict{Symbol,Any}(),
    withEsc::Bool = false,
) = mapreduce(
    ref -> filterreference(filt, ref..., props = props, withEsc = withEsc),
    append!,
    refs,
) 
