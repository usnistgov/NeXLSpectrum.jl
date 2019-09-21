using SparseArrays
using Polynomials
using LinearAlgebra

"""
    buildfilter(det::Detector, a::Float64=1.0, b::Float64=1.0)::AbstractMatrix{Float64}

Build a top-hat filter for the specified detector with the specified top and base parameters.
"""
function buildfilter(det::Detector, a::Float64=1.0, b::Float64=1.0, full::Bool=false)::AbstractMatrix{Float64}
    intersect(a,b,c,d) = max(0.0, min(b,d) - max(a,c)) # Length of intersection [a,b) with [c,d)
    filtint(e0, e1, minb, mina, maxa, maxb) =
        intersect(minb, mina, e0, e1)*(-0.5/(mina-minb)) +
        intersect(mina, maxa, e0, e1)/(maxa-mina) +
        intersect(maxa, maxb, e0, e1)*(-0.5/(maxb-maxa))
    #i,j,v=Array{Int}(),Array{Int}(),Array{Float64}()
    cc = channelcount(det)
    ans = zeros(Float64, (cc, cc))
    for ch1 in 1:cc
        ee = 0.5*(energy(ch1, det)+energy(ch1+1,det)) # midpoint of channel
        res = resolution(ee,det)
        ea = ( ee-0.5*a*res, ee+0.5*a*res)
        eb = ( ea[1]-0.5*b*res, ea[2]+0.5*b*res)
        if full
            chmin, chmax = max(1,channel(eb[1],det)), min(cc,channel(eb[2],det))
        else
            chmin, chmax = channel(eb[1],det), channel(eb[2],det)
        end
        if (chmin>=1) && (chmax<=cc) # fit the full filter
            sum=0.0 # Ensure the sum of the top-hat is precisely 0.0
            for ch2 in chmin:chmax-1
                fi = filtint(energy(ch2,det), energy(ch2+1,det), eb[1], ea[1], ea[2], eb[2])
                sum+=fi
                ans[ch1,ch2] = fi
            end
            ans[ch1,chmax]=-sum # Ensure the sum is zero...
        end
    end
    return sparse(ans)
end

"""
    FilteredDatum

Carries the minimal data necessary to support filter-fitting a single
region-of-interest (continguous range of channles) and computing
useful output statistics.
"""
struct FilteredDatum
    identifier # A way of identifying this filtered datum
    scale::Float64 # A dose or other scale correction factor
    roi::UnitRange{Int} # ROI for the raw data
    ffroi::UnitRange{Int} # ROI for the filtered data
    data::Vector{Float64} # Spectrum data over ffroi
    filtered::Vector{Float64} # Filtered spectrum data over ffroi
    covariance::AbstractMatrix{Float64} # Channel covariance
end

function Base.show(io::IO, fd::FilteredDatum)
    print(io, "FilteredDatum[",fd.identifier,"]")
    print(io," over ROC[",fd.roi.start,",",fd.roi.stop,"]")
end

"""
    FilterFit

A packet of data for fitting spectra.
"""
struct FilterFit
    filtered::Vector{FilteredDatum}
    filter::AbstractMatrix{Float64}
    FilterFit(filt::AbstractMatrix{Float64}) = new(Vector{FilteredDatum}(), filt)
end

function Base.show(io::IO, fp::FilterFit)
    println(io, "FilterFit with ",size(fp.filteredData)," items")
    for fd in fp.filtered
        println(io,fd)
    end
end

"""
    estimateBackground(data::AbstractArray{Float64}, channel::Int, width::Int=5, order::Int=2)

Returns the tangent to the a quadratic fit to the counts data centered at channel with width
"""
function estimateBackground(data::AbstractArray{Float64}, channel::Int, width::Int=5, order::Int=2)::Poly
    fit = polyfit(-width:width, data[max(1,channel-width):min(length(data),channel+width)],order)
    return Poly( [ fit(0), polyder(fit)(0) ] )
end

"""
    filter(id, spec::Spectrum, roi::UnitRange{Int},  filter::AbstractMatrix{Float64})::FilteredDatum

For filtering an ROI on a reference spectrum. Process a portion of a Spectrum with the specified filter.
"""
function Base.filter(id, spec::Spectrum, roi::UnitRange{Int},  filter::AbstractMatrix{Float64}, scale=1.0, tol::Float64 = 1.0e-4)::FilteredDatum
    data = counts(spec, Float64) # Extract the spectrum data
    # Determine tangents to the two background end points
    tangents = map(st->estimateBackground(data, st, 5, 2), (roi.start, roi.stop))
    # Replace the non-ROI channels with extensions of the tangent functions
    data[1:roi.start-1] = map(tangents[1], 1-roi.start:-1)
    data[roi.stop+1:end] = map(tangents[2], 1:length(data)-roi.stop)
    # Compute the filtered data
    filtered = filter*data
    maxval = maximum(filtered)
    roiff = findfirst(f->abs(f)>tol*maxval, filtered):findlast(f->abs(f)>tol*maxval, filtered)
    # max(d,1.0) is necessary to ensure the varianes are positive
    covar = filter*Diagonal(map(d->max(d,1.0),data))*transpose(filter)
    # Extract the range of non-zero filtered data
    trimmedcovar=sparse(covar[roiff,roiff])
    checkcovariance!(trimmedcovar)
    FilteredDatum(id, scale, roi, roiff, data[roiff], filtered[roiff], trimmedcovar)
end

"""
    filter(id, spec::Spectrum, filter::AbstractMatrix{Float64}, scale::Float64=1.0, tol::Float64 = 1.0e-4)::FilteredDatum

For filtering the unknown spectrum. Process the full Spectrum with the specified filter.
"""
function Base.filter(spec::Spectrum, filter::AbstractMatrix{Float64}, scale=1.0)::FilteredDatum
    data = counts(spec, Float64) # Extract the spectrum data
    # Compute the filtered data
    filtered = filter*data
    # max(d,1.0) is necessary to ensure the varianes are positive
    covar = filter*Diagonal(map(d->max(d,1.0),data))*transpose(filter)
    # Extract the range of non-zero filtered data
    roi = 1:findlast(f -> f ≠ 0.0,filtered)
    trimmedcovar=sparse(covar[roi,roi])
    checkcovariance!(trimmedcovar)
    FilteredDatum(UnknownLabel(spec), scale, roi, roi, data[roi], filtered[roi], trimmedcovar)
end


"""
    add!(ff::FilterFit, fd::FilteredDatum)

Add a FilteredDatum structure representing a range to be fit into
a FilterFit object representing the entire fit.
"""
add!(ff::FilterFit, fd::FilteredDatum) =
    push!(ff.filtered, fd)
"""
    extract(fd::FilteredDatum, roi::UnitRange{Int})

Extract the filtered data representing the specified range.
"""
function extract(fd::FilteredDatum, roi::UnitRange{Int})
    if (roi.start>fd.ffroi.stop) || (roi.stop < fd.ffroi.start)
        error("The FilteredDatum ",fd.identifier,ffroi," does not intersect "+roi)
    end
    # How many zeros do we need to add the beginning and end
    nz = max(1,1+fd.ffroi.start-roi.start):(length(roi)-max(0,roi.stop-fd.ffroi.stop))
    # Determine the range of channels to extract from fd.filter
    ffr = (1+max(0,roi.start-fd.ffroi.start)):min(length(fd.ffroi),1+roi.stop-fd.ffroi.start)
    # println("extract:: roi = ",roi,"  ff = ",fd.ffroi, "  nz = ",nz, "  ffr = ",ffr)
    data = zeros(Float64,length(roi))
    # Set the non-zero extent
    data[nz] = fd.filtered[ffr]
    return data
end

"""
    covariance(fd::FilteredDatum, roi::UnitRange{Int})

Like extract(fd,roi) except extracts the covariance matrix over the specified range of channels.
"""
function covariance(fd::FilteredDatum, roi::UnitRange{Int})
    if (roi.start>fd.ffroi.stop) || (roi.stop < fd.ffroi.start)
        error("The FilteredDatum ",fd.identifier,ffroi," does not intersect "+roi)
    end
    # How many zeros do we need to add the beginning and end
    nz = max(1,1+fd.ffroi.start-roi.start):(length(roi)-max(0,roi.stop-fd.ffroi.stop))
    # Determine the range of channels to extract from fd.filter
    ffr = (1+max(0,roi.start-fd.ffroi.start)):min(length(fd.ffroi),1+roi.stop-fd.ffroi.start)
    # println("covariance:: roi = ",roi,"  ff = ",fd.ffroi, "  nz = ",nz, "  ffr = ",ffr)
    cov = zeros(Float64,length(roi),length(roi))
    # Set the non-zero extent
    cov[nz,nz] = fd.covariance[ffr,ffr]
    return cov
end

"""
    ascontiguous(rois::AbstractArray{UnitRange{Int}})

Combine the array of ranges into a array in which intersecting ranges are
combined.  Returns a minimal set of contiguous ranges.
"""
function ascontiguous(rois::AbstractArray{UnitRange{Int}})
    join(roi1,roi2) = min(roi1.start,roi2.start):max(roi1.stop,roi2.stop)
    srois = sort(rois)
    res = [ srois[1] ]
    for roi in srois[2:end]
        if length(intersect(res[end],roi))>0
            res[end] = join(roi, res[end])
        else
            push!(res, roi)
        end
    end
    res
end

function fitcontiguousg(unk::FilteredDatum, ffs:: Array{FilteredDatum}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef,(length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:,i]=extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale/ffs[i].scale for i in eachindex(ffs)])
    return scale*glssvd(extract(unk,chs), A, covariance(unk, chs), xlbls)
end

function fitcontiguousp(unk::FilteredDatum, ffs:: Array{FilteredDatum}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef,(length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:,i]=extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale/ffs[i].scale for i in eachindex(ffs)])
    return scale*glspinv(extract(unk,chs), A, covariance(unk, chs), xlbls)
end

function fitcontiguousw(unk::FilteredDatum, ffs:: Array{FilteredDatum}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef,(length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:,i]=extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale/ffs[i].scale for i in eachindex(ffs)])
    return scale*wlssvd(extract(unk,chs), A, diag(covariance(unk,chs)), xlbls)
end

function fitcontiguouso(unk::FilteredDatum, ffs:: Array{FilteredDatum}, chs::UnitRange{Int})::UncertainValues
    # Build the fitting matrix
    A = Matrix{Float64}(undef,(length(chs), length(ffs)))
    for i in eachindex(ffs)
        A[:,i]=extract(ffs[i], chs)
    end
    # Build labels and scale
    xlbls = collect(ff.identifier for ff in ffs)
    scale = Diagonal([unk.scale/ffs[i].scale for i in eachindex(ffs)])
    return scale*olssvd(extract(unk,chs), A, 1.0, xlbls)
end

"""
    filterfit(unk::FilteredDatum, ffs::Array{FilteredDatum}, alg=fitcontiguousw)::UncertainValues

Filter fit the unknown against ffs, an array of FilteredDatum and return the result
as an UncertainValues object. Use generalized LLSQ fitting.
"""
function filterfit(unk::FilteredDatum, ffs::Array{FilteredDatum}, alg=fitcontiguousw)::UncertainValues
    println("Fitting: ",unk)
    join(roi1,roi2) = min(roi1.start,roi2.start):max(roi1.stop,roi2.stop)
    dup = copy(ffs)
    uvss=Array{UncertainValues,1}()
    while length(dup)>0
        ffs2=Array{FilteredDatum,1}()
        chs=1:1
        for ff in dup
            if length(ffs2) == 0
                push!(ffs2,ff)
                chs=ff.ffroi
            elseif length(intersect(ff.ffroi,chs)) > 0
                push!(ffs2,ff)
                chs = join(chs,ff.ffroi)
            end
        end
        println("References: ",ffs2)
        push!(uvss,alg(unk, ffs2, chs))
        for ff in ffs2
            splice!(dup,findfirst(d->d==ff,dup))
        end
    end
    println("\n")
    cat(uvss)
end

"""
    filteredresidual(fit::UncertainValues, unk::FilteredDatum, ffs::Array{FilteredDatum})::Vector{Float64}

Computes the difference between the best fit and the unknown filtered spectral data.
"""
filteredresidual(fit::UncertainValues, unk::FilteredDatum, ffs::Array{FilteredDatum})::Vector{Float64} =
     unk.filtered - mapreduce(ff->(value(ff.identifier,fit)*(ff.scale/unk.scale))*extract(ff,unk.ffroi), +, ffs)

function residual(fit::UncertainValues, unk::Spectrum, ffs::Array{Spectrum})::Vector{Float64}


end
