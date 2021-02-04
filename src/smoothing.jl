abstract type SpectrumFilter end

# From: https://gist.github.com/jiahao/b8b5ac328c18b7ae8a17


struct SavitzkyGolayFilter{M,N} <: SpectrumFilter end
@generated function (::SavitzkyGolayFilter{M,N})(data::AbstractVector{T}) where {M,N,T}
    # Create Jacobian matrix
    J = zeros(2M + 1, N + 1)
    for i = 1:2M+1, j = 1:N+1
        J[i, j] = (i - M - 1)^(j - 1)
    end
    e₁ = zeros(N + 1)
    e₁[1] = 1.0

    # Compute filter coefficients
    C = J' \ e₁
    # Evaluate filter on data matrix

    To = typeof(C[1] * one(T)) # Calculate type of output
    expr = quote
        n = size(data, 1)
        smoothed = zeros($To, n)
        @inbounds for i in eachindex(smoothed)
            smoothed[i] += $(C[M+1]) * data[i]
        end
        smoothed
    end

    for j = 1:M
        insert!(
            expr.args[6].args[3].args[2].args,
            1,
            :(
                if i - $j ≥ 1
                    smoothed[i] += $(C[M+1-j]) * data[i-$j]
                end
            ),
        )
        push!(
            expr.args[6].args[3].args[2].args,
            :(
                if i + $j ≤ n
                    smoothed[i] += $(C[M+1+j]) * data[i+$j]
                end
            ),
        )
    end
    return expr
end

"""
    apply(filt::SavitzkyGolayFilter, spec::Spectrum, applyLLD=false)

Applys a function to the channel data in `spec` (with/wo the low-level discriminator.)
The function can only be a function of the counts data and can not change the
energy scale or spectrum properties.  The result is a Spectrum.

Example:

    apply(SavitzkyGolayFilter{6,4}(),spec)
"""
function apply(filt::SavitzkyGolayFilter, spec::Spectrum, applyLLD = false)
    res = Spectrum(
        spec.energy,
        filt(counts(spec, eltype(spec), applyLLD)),
        copy(spec.properties),
    )
    res[:Name] = repr(filt)[1:end-1] * "$(spec[:Name]))"
    res[:Parent] = spec
    res
end
