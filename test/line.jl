using Test
using NeXLSpectrum
using Random


@testset "Line" begin
between(x, a, b) = x >= min(a, b) && x <= max(a, b)

N = 100
mt = MersenneTwister(0xBADF00C)
ci1s = CartesianIndex.(zip(rand(mt, -1000:1000, N), rand(mt, -1000:1000, N)))
ci2s = CartesianIndex.(zip(rand(mt, -1000:1000, N), rand(mt, -1000:1000, N)))

res = true
for (ci1, ci2) in zip(ci1s, ci2s)
    # println(ci1, " to ", ci2)
    drawline(ci1, ci2, false) do pt
        @assert between(pt[1], ci1[1], ci2[1]) "Point 1 failure $(pt[1]), $(ci1[1]), $(ci2[1])"
        res &= between(pt[1], ci1[1], ci2[1])
        @assert between(pt[2], ci1[2], ci2[2]) "Point 2 failure $(pt[2]), $(ci1[2]), $(ci2[2])"
        res &= between(pt[2], ci1[2], ci2[2])
        t1 = (pt[1] - ci1[1]) / (ci2[1] - ci1[1])
        t2 = (pt[2] - ci1[2]) / (ci2[2] - ci1[2])
        test1 = pt[1] - (ci1[1] + t2 * (ci2[1] - ci1[1]))
        test2 = pt[2] - (ci1[2] + t1 * (ci2[2] - ci1[2]))
        @assert abs(test1) <= 0.5 || abs(test2) <= 0.5 "Type 3 failure  $t1 $test1 $t2 $test2"
        res &= abs(test1) <= 0.5 || abs(test2) <= 0.5
    end
end
@test res

end
