"""
    olssvd(y::Vector{N}, a::AbstractMatrix{N}, w::AbstractMatrix{N}, xlabels::Vector{Label}, tol::N=convert(N,1.0e-10))::Vector{N}
Solves the ordinary least squares problem a⋅x = y for x using singular value decomposition for AbstractFloat-based types.
"""
function olssvd(y::AbstractVector{N}, a::AbstractMatrix{N}, sigma::N, xLabels::Vector{<:Label}, tol::N=convert(N,1.0e-10))::UncertainValues where N <: AbstractFloat
    f = svd(a)
    mins = tol*maximum(f.S)
    fs = [s > mins ? one(N)/s : zero(N) for s in f.S]
    # Note: w*a = f.U * Diagonal(f.S) * f.Vt
    genInv = f.V*Diagonal(fs)*transpose(f.U)
    #cov = sigma^2*f.V*(Diagonal(fs)*Diagonal(fs))*Transpose(f.V)
    cov = sigma^2*genInv*transpose(genInv)
    return uvs(xLabels, genInv*y, cov)
end

"""
    glspinv(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::Matrix{N}, xlabels::Vector{Label}, tol::N=convert(N,1.0e-10))::Vector{N}
Solves the generalized least squares problem a⋅x = y for x using singular value decomposition for AbstractFloat-based types.
"""
function glspinv(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::AbstractMatrix{N}, xLabels::Vector{<:Label}, tol::N=convert(N,1.0e-10))::UncertainValues where N <: AbstractFloat
    checkcovariance!(cov)
    ich = inv(cholesky(cov).L)
    genInv = pinv(ich*a, rtol=1.0e-6)
    return uvs(xLabels, genInv*(ich*y), genInv*transpose(genInv))
end

"""
    glssvd(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::Matrix{N}, xlabels::Vector{Label}, tol::N=convert(N,1.0e-10))::Vector{N}
Solves the generalized least squares problem a⋅x = y for x using singular value decomposition for AbstractFloat-based types.
"""
function glssvd(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::AbstractMatrix{N}, xLabels::Vector{<:Label}, tol::N=convert(N,1.0e-10))::UncertainValues where N <: AbstractFloat
    checkcovariance!(cov)
    ich = inv(cholesky(cov).L)
    olssvd(ich*y, ich*a, 1.0, xLabels, tol)
end


"""
    glssvd(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::AbstractVector{N}, xlabels::Vector{Label}, tol::N=convert(N,1.0e-10))::Vector{N}
Solves the generalized least squares problem a⋅x = y for x using singular value decomposition for AbstractFloat-based types.
"""
function wlssvd(y::AbstractVector{N}, a::AbstractMatrix{N}, cov::AbstractVector{N}, xLabels::Vector{<:Label}, tol::N=convert(N,1.0e-10))::UncertainValues where N <: AbstractFloat
    ich = Diagonal([sqrt(1.0/cv) for cv in cov])
    olssvd(ich*y, ich*a, 1.0, xLabels, tol)
end
