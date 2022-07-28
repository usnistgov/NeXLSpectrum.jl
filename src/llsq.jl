"""
    simple_linear_regression(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Tuple{Real, Real}

Peforms simple linear regression. Returns the `( slope, intercept, r, t )` of the unweighted best fit line through
the data `x` and `y`.  ( y = slope*x + intercept ), r is the correlation coefficient and t is the t-statistic.
"""
function simple_linear_regression(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    sx, sy, n = sum(x), sum(y), length(x)
    ssxy, ssxx, ssyy = (dot(x, y) - sx*sy/n), (dot(x, x) - sx*sx/n), (dot(y, y) - sy*sy/n)
    slope, r = ssxy / ssxx, ssxy/sqrt(ssxx*ssyy)
    return ( slope, (sy - slope*sx)/n, r, r/sqrt((1.0-r^2)/(n-2)) )
end
