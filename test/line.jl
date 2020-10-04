using Test
using NeXLSpectrum
using Random


between(x,a,b) = x >= min(a,b) && x<=max(a,b)

N = 100
mt = MersenneTwister(0xBADF00D)
ci1s = CartesianIndex.( zip(rand(mt, -1000:1000,N), rand(mt, -1000:1000,N)) )
ci2s = CartesianIndex.( zip(rand(mt, -1000:1000,N), rand(mt, -1000:1000,N)) )

for (ci1, ci2) in zip(ci1s, ci2s)
    # println(ci1, " to ", ci2)
    drawline( ci1, ci2, false) do pt
        @test between(pt[1], ci1[1], ci2[1])
        @test between(pt[2], ci1[2], ci2[2])
        t1 = (pt[1] - ci1[1])/(ci2[1]-ci1[1])
        t2 = (pt[2] - ci1[2])/(ci2[2]-ci1[2])
        @test abs(pt[1]  - (ci1[1] + t2*(ci2[1]-ci1[1]))) <= 0.5 || #
            abs(pt[2]  - (ci1[2] + t1*(ci2[2]-ci1[2]))) <= 0.5
    end
end
