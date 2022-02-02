using Test


@testset "Efficiency" begin
    sdd = SDDEfficiency(ModeledWindow(MoxtekAP33))
    sili = SiLiEfficiency(ModeledWindow(BerylliumWindow, thickness=10.0e-4))

    @test isapprox(efficiency(sdd, 300.0), 0.152, atol = 0.001)
    @test isapprox(efficiency(sdd, 3000.0), 0.753, atol = 0.001)
    @test isapprox(efficiency(sdd, 30000.0), 0.091, atol = 0.001)

    @test isapprox(efficiency(sili, 300.0), 0.000, atol = 0.001)
    @test isapprox(efficiency(sili, 3000.0), 0.957, atol = 0.001)
    @test isapprox(efficiency(sili, 30000.0), 0.485, atol = 0.001)
end
