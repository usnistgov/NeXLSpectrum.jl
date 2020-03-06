using Test
using NeXLCore
using NeXLSpectrum

@testset "Al2O3" begin
    dir = @__DIR__
    path = "$(dir)/K412 spectra/"
    al2o3 = readEMSA(path * "Al2O3 std.msa")
    @test isapprox(al2o3[:LiveTime],1172.19288,atol=0.00001)
    @test isapprox(al2o3[:BeamEnergy],20.0e3,atol=1.0)
    @test isapprox(al2o3[:ProbeCurrent],1.10989,atol=0.00001)
    @test isapprox(al2o3[:RealTime],1491.4828,atol=0.0001)
    @test isapprox(al2o3[:Elevation],deg2rad(35.0),atol=deg2rad(0.01))
    @test isapprox(material(al2o3[:Coating])[n"C"],1.0,atol=0.0001)
    @test isapprox(NeXLCore.thickness(al2o3[:Coating]),10.0*1.e-7,atol=1.0e-9)
    @test isapprox(al2o3[101],22290.0,atol=0.001)
    @test isapprox(energy(1, al2o3.energy),1.63032,atol=0.1)
    @test isapprox(energy(101, al2o3.energy),100.0*9.99856+1.63032,atol=0.1)
    @test isapprox(energy(1001, al2o3.energy),1000.0*9.99856+1.63032,atol=0.1)
    @test isapprox(energy(10001, al2o3.energy),10000.0*9.99856+1.63032,atol=0.1)
    @test channel(10000.0*9.99856+1.63032, al2o3.energy)==10001
    @test channel(1.63032-1.0, al2o3.energy)==0
    @test channel(100.0*9.99856+1.63032+1.0, al2o3.energy)==101
    @test channel(100.0*9.99856+1.63032-1.0, al2o3.energy)==100
    @test al2o3[:Name]=="Al2O3 std"
    @test isapprox(al2o3[:Composition][n"Al"],0.529251,atol=0.00001)
    @test isapprox(al2o3[:Composition][n"O"],0.470749,atol=0.00001)
    @test isapprox(al2o3[:WorkingDistance],1.7,atol=0.0001)
    @test al2o3[:Owner]=="Unknown"
end
