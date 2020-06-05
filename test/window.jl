using Test

@testset "Windows" begin
    ap33 = AP33Model()
    @test repr(ap33)=="Moxtek AP3.3 model"
    ap5 = AP5Model()
    @test repr(ap5)=="Moxtek AP5 model"
    ap33t = AP33Tabulation()
    @test repr(ap33t)=="Moxtek AP3.3"
    ap5t = AP5Tabulation()
    @test repr(ap5t)=="Moxtek AP5"
    be = Beryllium()
    @test repr(be)=="5.0 μm Be window"
    ac1 = AmptekC1()
    @test repr(ac1)=="Amptek C1"
    ac2 = AmptekC2()
    @test repr(ac2)=="Amptek C2"
    now = NoWindow() # 100% transmission
    @test repr(now)=="No window"

    @test transmission(now, 300.0) == 1.0
    @test isapprox(transmission(ap33, 300.0),0.197,atol=0.001)
    @test isapprox(transmission(ap33t, 300.0),0.342,atol=0.001)
    @test isapprox(transmission(ap5t, 300.0),0.132,atol=0.001)
    @test isapprox(transmission(ap5, 300.0),0.200,atol=0.001)
    @test isapprox(transmission(be, 300.0), 3.6e-7, atol=0.1e-7)
    @test isapprox(transmission(ac1, 300.0), 0.0702, atol=0.001)
    @test isapprox(transmission(ac2, 300.0), 0.5198, atol=0.001)


    @test transmission(now, 3000.0) == 1.0
    @test isapprox(transmission(ap33, 3000.0),0.760,atol=0.001)
    @test isapprox(transmission(ap33t, 3000.0),0.772,atol=0.001)
    @test isapprox(transmission(ap5t, 3000.0),0.772,atol=0.001)
    @test isapprox(transmission(ap5, 3000.0),0.773,atol=0.001)
    @test isapprox(transmission(be, 3000.0), 0.982, atol=0.001)
    @test isapprox(transmission(ac1, 3000.0), 0.865, atol=0.001)
    @test isapprox(transmission(ac2, 3000.0), 0.928, atol=0.001)

    @test isapprox(transmission(ap33, 30000.0),0.978,atol=0.001)
    @test isapprox(transmission(ap33t, 30000.0),0.993,atol=0.001)
    @test isapprox(transmission(ap5t, 30000.0),0.999,atol=0.001)
    @test isapprox(transmission(ap5, 30000.0),0.999,atol=0.001)
    @test isapprox(transmission(be, 30000.0), 0.9999, atol=0.001)
    @test isapprox(transmission(ac1, 30000.0), 0.9999, atol=0.0001)
    @test isapprox(transmission(ac2, 30000.0), 0.9999, atol=0.0001)
end
