using Test

@testset "Windows" begin
    ap33 = ModeledWindow(MoxtekAP33())
    @test repr(ap33) == "Moxtek AP3.3 - Modeled"
    ap5 = ModeledWindow(MoxtekAP5())
    @test repr(ap5) == "Moxtek AP5 - Modeled"
    ap33t = TabulatedWindow(MoxtekAP33())
    @test repr(ap33t) == "Moxtek AP3.3 - Tabulated"
    ap5t = TabulatedWindow(MoxtekAP5())
    @test repr(ap5t) == "Moxtek AP5 - Tabulated"
    be = ModeledWindow(BerylliumWindow(5.0e-4))
    @test repr(be) == "5.0 μm Beryllium - Modeled"
    ac1 = ModeledWindow(AmetekC1())
    @test repr(ac1) == "AMETEK C1 Si₃N₄ - Modeled"
    ac2 = ModeledWindow(AmetekC2())
    @test repr(ac2) == "AMETEK C2 Si₃N₄ - Modeled"
    now = NoWindow() # 100% transmission
    @test repr(now) == "No window"

    @test transmission(now, 300.0) == 1.0
    @test isapprox(transmission(ap33, 300.0), 0.197, atol = 0.001)
    @test isapprox(transmission(ap33t, 300.0), 0.342, atol = 0.001)
    @test isapprox(transmission(ap5t, 300.0), 0.132, atol = 0.001)
    @test isapprox(transmission(ap5, 300.0), 0.200, atol = 0.001)
    @test isapprox(transmission(be, 300.0), 3.6e-7, atol = 0.1e-7)
    @test isapprox(transmission(ac1, 300.0), 0.0702, atol = 0.001)
    @test isapprox(transmission(ac2, 300.0), 0.5198, atol = 0.001)


    @test transmission(now, 3000.0) == 1.0
    @test isapprox(transmission(ap33, 3000.0), 0.760, atol = 0.001)
    @test isapprox(transmission(ap33t, 3000.0), 0.772, atol = 0.001)
    @test isapprox(transmission(ap5t, 3000.0), 0.772, atol = 0.001)
    @test isapprox(transmission(ap5, 3000.0), 0.773, atol = 0.001)
    @test isapprox(transmission(be, 3000.0), 0.982, atol = 0.001)
    @test isapprox(transmission(ac1, 3000.0), 0.865, atol = 0.001)
    @test isapprox(transmission(ac2, 3000.0), 0.928, atol = 0.001)

    @test isapprox(transmission(ap33, 30000.0), 0.978, atol = 0.001)
    @test isapprox(transmission(ap33t, 30000.0), 0.993, atol = 0.001)
    @test isapprox(transmission(ap5t, 30000.0), 0.999, atol = 0.001)
    @test isapprox(transmission(ap5, 30000.0), 0.999, atol = 0.001)
    @test isapprox(transmission(be, 30000.0), 0.9999, atol = 0.001)
    @test isapprox(transmission(ac1, 30000.0), 0.9999, atol = 0.0001)
    @test isapprox(transmission(ac2, 30000.0), 0.9999, atol = 0.0001)
end

