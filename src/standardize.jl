using NeXLUncertainties

"""
    StandardLabel

A ReferenceLabel that represents a reference spectrum or reference properties associated with a set of 
characteristic x-rays (CharXRay) objects over a contiguous range of spectrum channels that is to be
used as a standard.
"""
struct StandardLabel <: ReferenceLabel
    spectrum::SpectrumOrProperties # The spectrum used as the reference for fitting...
    # roi::UnitRange{Int} # Not available in k-ratios
    xrays::Vector{CharXRay}
    standard::Material
    hash::UInt
    function StandardLabel(cxrl::CharXRayLabel, stdMat::Material)
        new(cxrl.spectrum, cxrl.xrays, stdMat, hash(cxrl.spectrum, hash(cxrl.xrays, hash(stdMat, UInt(0xBADF00D5)))))
    end
    function StandardLabel(props::Dict{Symbol,Any}, xrays::Vector{CharXRay}, stdMat::Material)
        new(props, xrays, stdMat, hash(props, hash(xrays, hash(stdMat, UInt(0xBADF00D5)))))
    end

end

"""
    matches(cxrl::CharXRayLabel, std::KRatio)
    matches(cxrl::CharXRayLabel, std::StandardLabel)

Does the measured `CharXRayLabel` match the standard?
"""
function matches(cxrl::CharXRayLabel, std::KRatio)::Bool
    cp, sp = NeXLCore.properties(cxrl), std.stdProps
    return issetequal(cxrl.xrays, std.xrays) && # Same X-rays
        isequal(cp[:Composition], sp[:Composition]) && # Same reference material
        isapprox(cp[:BeamEnergy], sp[:BeamEnergy], rtol=0.001) && # Same beam energy
        isapprox(cp[:TakeOffAngle], sp[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off-angle
end
function matches(cxrl::CharXRayLabel, std::StandardLabel)::Bool
    cp, sp = NeXLCore.properties(cxrl), NeXLCore.properties(std)
    return issetequal(cxrl.xrays, std.xrays) && # Same X-rays
        isequal(cp[:Composition], sp[:Composition]) && # Same reference material
        isapprox(cp[:BeamEnergy], sp[:BeamEnergy], rtol=0.001) && # Same beam energy
        isapprox(cp[:TakeOffAngle], sp[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off-angle
end

function Base.show(io::IO, cl::StandardLabel) 
    comp = composition(cl)
    compname = isnothing(comp) ? "Unspecified" : name(comp)
    print(io,"k_std[$(name(cl.xrays)), $(name(cl.standard)) fit by $compname]")
end

"""
`StandardizeModel` is designed to combine k-ratio measurements relative to common reference spectra
to re-standardize a microanalysis measurement.  As a `MeasurementModel` it combines k-ratio data
while maintaining covariance relationships between the measured values.
"""
struct StandardizeModel <: MeasurementModel 
    measured::Vector{CharXRayLabel}
    standards::Vector{StandardLabel}

    StandardizeModel(measured::AbstractVector{CharXRayLabel}, standards::AbstractVector{StandardLabel}) = new(measured, standards)
end

function NeXLUncertainties.compute(
    sm::StandardizeModel,
    inputs::LabeledValues,
    withJac::Bool,
)::MMResult
    outputs, results = Array{CharXRayLabel}(undef, length(sm.measured)), zeros(Float64, length(sm.measured))
    for (i, meas) in enumerate(sm.measured)
        j = findfirst(std->matches(meas, std), sm.standards)
        if (!isnothing(j)) && (inputs[sm.standards[j]] > 0.0) # Re-standardized
            std = sm.standards[j] 
            sp = copy(NeXLCore.properties(std))
            sp[:Composition] = std.standard
            outputs[i] = CharXRayLabel(sp, meas.roi, meas.xrays) # restandardized label
            results[i] = inputs[meas] / inputs[std] 
        else
            outputs[i] = meas
            results[i] = inputs[meas]
        end
    end
    if withJac
        jac = zeros(Float64, length(sm.measured), length(sm.measured)+length(sm.standards))
        for (i, meas) in enumerate(sm.measured)
            j = findfirst(std->matches(meas, std), sm.standards)
            if (!isnothing(j)) &&  (inputs[sm.standards[j]] > 0.0) # Re-standardized
                std = sm.standards[j]
                jac[i, indexin(inputs, meas)] = 1.0 / inputs[std]
                jac[i, indexin(inputs, std)] = -inputs[meas]/(inputs[std]^2)
            else
                jac[i, indexin(inputs, meas)] = 1.0
            end
        end
    end
    return (LabeledValues(outputs, results), jac)
end

"""
    NeXLCore.standardize(ffr::FilterFitResult{T}, standard::FilterFitResult{T}, material::Material, els=elms(material))::FilterFitResult{T}
    NeXLCore.standardize(ffrs::Vector{FilterFitResult{T}}, standard::FilterFitResult{T}, material::Material, els=elms(material))::Vector{FilterFitResult{T}}
    NeXLCore.standardize(ffr::FilterFitResult{T}, standards::AbstractArray{KRatio})::FilterFitResult{T}
    NeXLCore.standardize(ffrs::Vector{FilterFitResult{T}}, standards::AbstractArray{KRatio})::Vector{FilterFitResult{T}}

Apply the standard `KRatio`s to the `FilterFitResult` producing a re-standardized `FilterFitResult`.
"""
function NeXLCore.standardize(ffr::FilterFitResult{T}, standard::FilterFitResult{T}, material::Material, els=elms(material)) where { T<: AbstractFloat }
    stds = filter(labels(standard.kratios)) do lbl
        (lbl isa CharXRayLabel) && (element(lbl) in els)
    end
    ext = extract(standard.kratios, stds)
    rext = uvs([ StandardLabel(l, material) for l in labels(ext)], ext.values, ext.covariance)
    # Problem here with equivalent CharXRayLabel(s) for unknown and standard...
    inp = cat(ffr.kratios, rext)
    stdize = StandardizeModel(Vector{CharXRayLabel}(labels(ffr.kratios)), Vector{StandardLabel}(labels(rext)))
    ptob = Deferred() do
        pbs, res = ffr.peakback, Dict{ReferenceLabel,NTuple{3,T}}()
        for (meas, pb) in pbs()
            i = findfirst(std->matches(meas, std), sm.standards)
            if !isnothing(i)
                std = sm.standards[i]
                sp = copy(NeXLCore.properties(std))
                sp[:Composition] = std.standard
                res[CharXRayLabel(sp, meas.roi, meas.xrays)] = pb
            else
                res[meas] = pb
            end
        end
        return res
    end
    return FilterFitResult{T}(
            ffr.label, # Same
            stdize(inp), # remapped
            ffr.roi, # Same
            ffr.raw, # Same
            ffr.residual, # Same - residual is due to references
            ptob 
        )
end

function NeXLCore.standardize(ffrs::AbstractVector{FilterFitResult{T}}, standard::FilterFitResult{T}, material::Material, els=elms(material)) where { T<: AbstractFloat }
    map(ffrs) do ffr
        standardize(ffr, standard, material, els)
    end
end

function NeXLCore.standardize(ffr::FilterFitResult{T}, standards::AbstractArray{KRatio}) where { T<: AbstractFloat }
    @assert all(std->isstandard(std), standards) "One or more required property is missing from a standard (for details, see the above warnings.)" 
    stdis = map(cxrl->findfirst(std->matches(cxrl, std), standards), labels(ffr.kratios))
    stds = uvs(map(filter(i->!isnothing(i), stdis)) do i
        std = standards[i]
        StandardLabel(std.stdProps, std.xrays, std.unkProps[:Composition])=>uv(std.kratio)
    end...)
    return if length(stds)>0
        stdize = StandardizeModel(Vector{CharXRayLabel}(labels(ffr.kratios)), Vector{StandardLabel}(labels(stds)))
        inp = cat(ffr.kratios, stds)
        ptob = Deferred() do
            pbs, res = ffr.peakback, Dict{ReferenceLabel,NTuple{3,T}}()
            for (meas, pb) in pbs()
                i = findfirst(std->matches(meas, std), sm.standards)
                if !isnothing(i)
                    std = sm.standards[i]
                    sp = copy(NeXLCore.properties(std))
                    sp[:Composition] = std.standard
                    res[CharXRayLabel(sp, meas.roi, meas.xrays)] = pb
                else
                    res[meas] = pb
                end
            end
            return res
        end
        FilterFitResult{T}(
            ffr.label,
            stdize(inp),
            ffr.roi,
            ffr.raw,
            ffr.residual,
            ptob
        )
    else
        @warn "No suitable standards were provided to standardize this measurement."
        ffr
    end
end

function NeXLCore.standardize(ffrs::AbstractVector{FilterFitResult{T}}, standards::AbstractArray{KRatio}) where { T<: AbstractFloat }
    map(ffrs) do ffr
        standardize(ffr, standards)
    end
end