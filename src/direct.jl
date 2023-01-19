# Implements a "direct fit" of reference spectra to an unknown spectrum
# This involves background modeling each spectrum - unknown and references - to remove the 
# signal from the continuum.  Then the remaining signal is assumed to be characteristic only
# so we can perform a direct linear least squares fit to determine the k-ratios.

"""
A `DirectReference` represents the data extracted from a reference spectrum necessary to perform a 
"direct" (background-eliminated direct fit)."""
struct DirectReference
    label::CharXRayLabel
    scale::Float64
    roi::UnitRange{Int}
    data::Vector{Float64}
end

Base.show(io::IO, dr::DirectReference) = print(io, "DirectReference($(dr.label))")

"""
A `DirectReferences` is a packet of `DirectReference` and detector information.  It can be used to
fit the corresponding peaks in an "unknown" spectrum.
"""
struct DirectReferences
    references::Vector{DirectReference}
    detector::Detector
    response::Matrix{Float64}
end

struct DirectRefInit
    element::Element
    spectrum::Spectrum
    material::Material
end

direct(elm::Element, spec::Spectrum, mat::Material) = DirectRefInit(elm, spec, mat)

"""
    references(refs::Vector{DirectRefInit}, det::Detector, resp::Matrix{Float64}; minE=0.5e3 )

Use along with `direct(...)` to build a collection of direct fitting references.

Example:

    > detu = matching(unk, 132.0, 110)
    > resp = detectorresponse(detu, SDDEfficiency(ModeledWindow(MoxtekAP33())))

    > drefs = references([ 
          direct(n"O", stds[1], mat"Al2O3"),
          direct(n"Al", stds[1], mat"Al2O3"),
          direct(n"Ba", stds[2], mat"BaF2"),
          direct(n"Ca", stds[3], mat"CaF2"),
          direct(n"Fe", stds[4], mat"Fe"),
          direct(n"Si", stds[5], mat"Si") ],
          detu, resp
      )
"""
function references(refs::Vector{DirectRefInit}, det::Detector, resp::Matrix{Float64}; minE=0.5e3 )
    function direct(elm::Element, spec::Spectrum, mat::Material)
        allElms = collect(elms(mat))
        @assert elm in allElms "$elm is not contained in $mat."
        @assert haskey(spec, :BeamEnergy) "The spectrum must define the :BeamEnergy property."
        lbls = NeXLSpectrum.charXRayLabels(spec, elm, allElms, det, spec[:BeamEnergy])
        map(lbls) do lbl
            charOnly = counts(lbl.spectrum, lbl.roi) - counts(fittedcontinuum(lbl.spectrum, det, resp, mode=:Local, minE=minE), lbl.roi)
            DirectReference(lbl, dose(lbl.spectrum), lbl.roi, charOnly)
        end
    end
    drefs = mapreduce(append!, refs) do ref
        direct(ref.element, ref.spectrum, ref.material)
    end
    return DirectReferences(drefs, det, resp)
end


Base.show(io::IO, drs::DirectReferences) = print(io, "DirectReferences(\n"*join( ("\t$(dr.label)" for dr in drs.references), ",\n")*"\n)")

"""
`DirectFitResult` contains the result of a direct fit of a `DirectReferences` to an unknown spectrum.
"""
struct DirectFitResult{T <: AbstractFloat} <: FitResult
    label::UnknownLabel
    kratios::UncertainValues
    residual::Spectrum{T}
    continuum::Spectrum{T}
    peakback::Dict{ReferenceLabel,NTuple{3,T}}
end

Base.show(io::IO, dfr::DirectFitResult) = print(io, "DirectFitResult($(dfr.label))")
residual(dfr::DirectFitResult) = dfr.residual
continuum(dfr::DirectFitResult) = dfr.continuum

"""
    fit_spectrum(unk::Spectrum{T}, drefs::DirectReferences)::DirectFitResult{T} where { T <: Real }

Fits the direct references to the unknown spectrum returning a `DirectFitResult` struct containing
k-ratios and other output quantities from the fit.

Surprisingly, this function requires that `unk[:Composition]` is defined as a `Material` with an
estimate of the composition of the material from which `unk` was collected.  This information is
necessary to model the continuum background.  Yes, this seems circular (and it is.)  One option
is to perform a filter-fit first, quantify the filter fit and then use this as your estimate.
"""
function fit_spectrum(unk::Spectrum{T}, drefs::DirectReferences) where { T <: Real }
    @assert haskey(unk.properties, :Composition) "Please define the :Composition property for `unk` with a `Material`."
    # Create blocks of contiguous peaks to fit in a group
    function contiguous(srefs)
        res = Vector{DirectReference}[ ]
        for (i, sref) in enumerate(sort(srefs.references, lt=(a,b)->a.roi.start < b.roi.start))
            if (i==1) || (sref.roi.start > res[end][end].roi.stop)
                push!(res, [ sref ])  # start new grouping
            else
                push!(res[end], sref) # append to last grouping
            end
        end
        return res
    end
    continuum = fittedcontinuum(unk, drefs.detector, drefs.response, mode=:Local, minE=0.5e3)
    charonly = unk - continuum
    du = dose(unk)
    # Perform the fit (returns `UncertainValues`)
    ks = cat(map(contiguous(drefs)) do block
        # Determine the fit ROI
        full = minimum(dr.roi.start for dr in block):maximum(dr.roi.stop for dr in block)
        # Extract these channels from charonly and unk and scale to 1 nA⋅s
        y, cov = counts(charonly, full) / du, counts(unk, full) / (du^2)
        # Build the model matrix
        a = zeros(Float64, length(y), length(block))
        for (i,dr) in enumerate(block)
            a[dr.roi.start-full.start+1:dr.roi.stop-full.start+1, i] = dr.data / dr.scale
        end
        # Perform a weighed least squares fit using the pseudo-inverse
        w = Diagonal([sqrt(1.0 / cv) for cv in cov])
        genInv = pinv(w * a, rtol = 1.0e-6)
        return uvs(getproperty.(block, :label), genInv * w * y, genInv * transpose(genInv))
    end)
    # Compute the residual spectrum
    function residual(unk, drefs, ks)
        res, du, rprops = copy(unk.counts), dose(unk), copy(unk.properties)
        for dr in drefs.references
            res[dr.roi] -= (value(ks, dr.label) * du/dr.scale) * dr.data
        end
        rprops[:Name] = "Residual[$(unk[:Name])]"
        return Spectrum(unk.energy, res, rprops)
    end
    DirectFitResult(
        UnknownLabel(unk),
        ks,
        residual(unk, drefs, ks),
        continuum,
        Dict{ReferenceLabel,NTuple{3,T}}()
    )
end