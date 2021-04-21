using NeXLUncertainties

"""
`StandardizeModel` is designed to combine k-ratio measurements relative to common reference spectra
to re-standardize a microanalysis measurement.  As a `MeasurementModel` it combines k-ratio data
while maintaining covariance relationships between the measured values.
"""
struct StandardizeModel <: MeasurementModel 
    measured::Vector{ReferenceLabel}
    standards::Vector{CharXRayLabel}

    StandardizeModel(measured::AbstractVector{ReferenceLabel}, standards::AbstractVector{CharXRayLabel}) = new(measured, standards)
end


function standardFor(meas::Label, sm::StandardizeModel)
    findfirst(sm.standards) do std
        (meas.roi===std.roi) && issetequal(meas.xrays, std.xrays)
    end
end

function NeXLUncertainties.compute(
    sm::StandardizeModel,
    inputs::LabeledValues,
    withJac::Bool,
)::MMResult
    outputs, results = Array{CharXRayLabel}(undef, length(sm.measured)), zeros(Float64, length(sm.measured))
    for (i, meas) in enumerate(sm.measured)
        j = standardFor(meas, sm)
        if (!isnothing(j)) && (inputs[sm.standards[j]] > 0.0) # Re-standardized
            std = sm.standards[j] 
            sp = copy(properties(std))
            @assert haskey(sp, :Composition) && haskey(sp, :BeamEnergy) && haskey(sp, :TakeOffAngle) 
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
            j = standardFor(meas, sm)
            if (!isnothing(j)) &&  (inputs[sm.standards[j]] > 0.0) # Re-standardized
                std = sm.standards[j]
                jac[i, indexin(meas, inputs)] = 1.0 / inputs[std]
                jac[i, indexin(std, inputs)] = -inputs[meas]/(inputs[std]^2)
            else
                jac[i, indexin(meas, inputs)] = 1.0
            end
        end
    end
    return (LabeledValues(outputs, results), jac)
end

"""
    matches(cxrl::CharXRayLabel, std::KRatio)
    matches(cxrl::CharXRayLabel, std::CharXRayLabel)

Does the measured `CharXRayLabel` match the standard?
"""
function matches(cxrl::CharXRayLabel, std::KRatio)::Bool
    cp, sp = properties(cxrl), std.stdProps
    return issetequal(cxrl.xrays, std.lines) && # Same X-rays
        isequal(cp[:Composition], sp[:Composition]) && # Same reference material
        isapprox(cp[:BeamEnergy], sp[:BeamEnergy], rtol=0.001) && # Same beam energy
        isapprox(cp[:TakeOffAngle], sp[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off-angle
end
function matches(cxrl::CharXRayLabel, std::CharXRayLabel)::Bool
    cp, sp = properties(cxrl), properties(std)
    return issetequal(cxrl.xrays, std.xrays) && # Same X-rays
        isequal(cp[:Composition], sp[:Composition]) && # Same reference material
        isapprox(cp[:BeamEnergy], sp[:BeamEnergy], rtol=0.001) && # Same beam energy
        isapprox(cp[:TakeOffAngle], sp[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off-angle
end

function __remap_peaktoback(pbs::Dict{<:ReferenceLabel,NTuple{2,Float64}}, sm::StandardizeModel)::Dict{ReferenceLabel,NTuple{2,Float64}}
    res = Dict{ReferenceLabel,NTuple{2,Float64}}()
    for (meas, pb) in pbs
        i = standardFor(meas, sm)
        lbl = isnothing(i) ? meas : CharXRayLabel(copy(properties(sm.standards[i])), meas.roi, meas.xrays)
        res[lbl] = pb
    end
    return res
end

"""
    NeXLCore.standardize(ffr::FilterFitResult, standard::FilterFitResult, material::Material, els=elms(material))::FilterFitResult
    NeXLCore.standardize(ffr::FilterFitResult, standards::Pair{<:Material,FilterFitResult}...)::FilterFitResult
    NeXLCore.standardize(ffr::FilterFitResult, standards::AbstractArray{KRatio})::FilterFitResult

Apply the standard `KRatio`s to the `FilterFitResult` producing a re-standardized `FilterFitResult`.
"""
function NeXLCore.standardize(ffr::FilterFitResult, standard::FilterFitResult, material::Material, els=elms(material))::FilterFitResult
    stds = filter(labels(standard.kratios)) do lbl
        properties(lbl)[:Composition] = material
        (lbl isa CharXRayLabel) && (element(lbl) in els)
    end
    stdize = StandardizeModel(convert(Vector{ReferenceLabel},labels(ffr.kratios)), convert(Vector{CharXRayLabel},stds))
    return FilterFitResult(
            ffr.label, # Same
            stdize(cat(ffr.kratios, extract(stds, standard.kratios))), # remapped
            ffr.roi, # Same
            ffr.raw, # Same
            ffr.residual, # Same - residual is due to references
            __remap_peaktoback(ffr.peakback, stdize) # Remapped to standard labels
        )
end

function NeXLCore.standardize(ffr::FilterFitResult, standards::Pair{<:Material,FilterFitResult}...)::FilterFitResult
    res = ffr
    for (std_mat, std_ffr) in standards
        res = standardize(res, std_ffr, std_mat)
    end
    return res
end

function NeXLCore.standardize(ffr::FilterFitResult, standards::AbstractArray{KRatio})::FilterFitResult
    @assert all(std->suitable_as_standard(std), standards) "One or more required property is missing from a standard (for details, see the above warnings.)" 
    stdis = map(cxrl->findfirst(std->matches(cxrl, std), standards), labels(ffr.kratios))
    stds = map(filter(i->!isnothing(i), stdis)) do i
        CharXRayLabel(std.unkProps, cxrl.roi, cxrl.xrays)=>uv(std.kratio)
    end
    return if length(stds)>0
        meas = labels(ffr.kratios)
        stdize = StandardizeModel(meas, stds)
        inp = cat(ffr.kratios, uvs( stds... ))
        FilterFitResult(
            ffr.label,
            stdize(inp),
            ffr.roi,
            ffr.raw,
            ffr.residual,
            __remap_peaktoback(ffr.peakback, sm)
        )
    else
        @warn "No suitable standards were provided to standardize this measurement."
        ffr
    end
end


