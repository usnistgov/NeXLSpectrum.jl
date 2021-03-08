using Test
using NeXLSpectrum
using NeXLCore
using Statistics
using Distributions

@testset "Spectrum" begin
    les = LinearEnergyScale(-495.0, 5.0)

    @testset "Energy Scale" begin
        @test NeXLSpectrum.energy(200, les) == -495.0 + 5.0 * (200 - 1)
        @test NeXLSpectrum.energy(1, les) == -495.0 + 5.0 * (1 - 1)
        @test NeXLSpectrum.energy(500, les) == -495.0 + 5.0 * (500 - 1)

        @test channel(NeXLSpectrum.energy(200, les) + 1.0, les) == 200
        @test channel(NeXLSpectrum.energy(200, les) - 1.0, les) == 199
        @test channel(NeXLSpectrum.energy(200, les), les) == 200

        @test channel(600.0, les) == 220
        @test channel(1600.0, les) == round((1600.0 + 495.0) / 5.0) + 1
    end

    mnk = MnKaResolution(130.0)

    @testset "Resolution" begin
        @test resolution(5898.7, mnk) == 130.0
        @test isapprox(resolution(500.0, mnk), 60.0, atol = 1.0)
    end

    @testset "Spectrum" begin
        s1 = Spectrum(
            LinearEnergyScale(0.0, 10.0),
            collect(1:2048),
            Dict{Symbol,Any}([(:BeamEnergy, 20.0e3), (:LiveTime, 60.0)]),
        )

        @test s1[:LiveTime] == 60.0
        @test s1[100] == 100
        @test s1[1000] == 1000

        @test s1[100:10:150] == [100, 110, 120, 130, 140, 150]
        s1[100] = 200
        s1[101:103] = [300, 400, 500]

        @test s1[99] == 99
        @test s1[100] == 200
        @test s1[101] == 300
        @test s1[102] == 400
        @test s1[103] == 500
        @test s1[104] == 104

        s1[:RealTime] = 90.0
        @test s1[:RealTime] == 90.0
    end

    @testset "SimpleEDS" begin
        eds = simpleEDS(4096, 5.0, -495.0, 125.0)
        @test channelcount(eds) == 4096
        @test isapprox(resolution(10000.0, eds), 160.0, atol = 1.0)
        @test NeXLSpectrum.energy(2000, eds) == NeXLSpectrum.energy(2000, les)
        @test channel(NeXLSpectrum.energy(2000, eds), eds) == 2000
    end
    @testset "BasicEDS" begin
        det3 = BasicEDS(
            4095,
            0.0,
            10.0,
            126.0,
            10,
            Dict(
                Shell(1) => n"Be",
                Shell(2) => n"Sc",
                Shell(3) => n"Cs",
                Shell(4) => n"Pu",
            ),
        )
        @test all(ch -> channel(energy(ch, det3), det3) == ch, 1:4096)
        @test resolution(energy(n"Mn K-L3"), det3) == 126.0
        @test isapprox(resolution(energy(n"C K-L2"), det3), 45.867, atol = 0.001)
        @test isvisible(n"C K-L2", det3)
        @test !isvisible(n"Ca L3-M3", det3)
        @test length(isvisible(characteristic(n"Sc", ltransitions), det3)) >= 7
        @test !isvisible(n"Xe M3-N5", det3)
        @test isvisible(n"Cs M3-N5", det3)
        @test length(isvisible(characteristic(n"U", ntransitions), det3)) == 0
        @test length(isvisible(characteristic(n"Pu", ntransitions), det3)) == 0 # Database doesn't contain any...
        @test channel(n"Mn K-L3", det3) == channel(energy(n"Mn K-L3"), det3)
    end
    @testset "readEMSA" begin
        sp = loadspectrum(joinpath(@__DIR__, "K412 spectra", "Al2O3 std.msa"))
        @test sp[:BeamEnergy] == 20.0e3
        @test sp[:ProbeCurrent] == 1.10989
        @test sp[:Elevation] == deg2rad(35.0)
        @test sp[:LiveTime] == 1172.19288
        @test sp[:RealTime] == 1491.4828
        @test sp[:Owner] == "Unknown"
        @test sp[:Name] == "Al2O3 std"
        comp = sp[:Composition]

        @test isapprox(comp[n"O"], 0.470749, atol = 1.0e-6)
        @test isapprox(comp[n"Al"], 0.529251, atol = 1.0e-6)
        @test comp[n"Zr"] == 0.0
        @test density(comp) == 4.00

        chs = channel(1335.0, sp):channel(1688.0, sp)
        @test integrate(sp, chs) == 33143490

        low, high = 120:133, 172:180
        @test isapprox(
            NeXLUncertainties.value(integrate(sp, low, chs, high)),
            3.22027e7,
            rtol = 3.0e-3,
        )

        ss = [subsample(sp, 0.1) for _ = 1:10]

        @test isapprox(
            sp[:LiveTime],
            sum(ss[i][:LiveTime] for i in eachindex(ss)),
            rtol = 1.0e-9,
        )
        @test isapprox(
            sp[:RealTime],
            sum(ss[i][:RealTime] for i in eachindex(ss)),
            rtol = 1.0e-9,
        )

        ch = channel(1450.0, sp)
        @test isapprox(
            sum(ss[i][ch] for i in eachindex(ss)),
            sp[ch],
            atol = 4.0 * sqrt(sp[ch]),
        ) # can fail occasionally...

    end
    @testset "SLR" begin
        x = [ 5.7, 2.5, 3.8, 2.8, 4.6, 1.6, 1.6, 1.4 ]
        y = [ 18, 19, 20, 19, 17, 32, 29, 30 ]
        slr = NeXLSpectrum.simple_linear_regression(x, y)
        @test isapprox(slr[1],-3.181, atol=1.0e-3)
        @test isapprox(slr[2], 32.542, atol=1.0e-3)
    end
    @testset "Duane-Hunt" begin
        fns = ( "Al2O3 std", "CaF2 std", "Fe std", "MgO std", "SiO2 std", ("III-E K412[$i][4]" for i in 0:4)...)
        specs = map(fn->loadspectrum(joinpath(@__DIR__,"K412 spectra",fn*".msa")),fns)
        @test all(isapprox.(NeXLSpectrum.duane_hunt.(specs), 20.0e3, atol=100.0))
    end

    @testset "sigma" begin
        specs = loadspectrum.(joinpath(@__DIR__,"K412 spectra","III-E K412[$i][4].msa") for i in 0:4)
        ss = [ sigma(specs[i], specs, 1:2000) for i in eachindex(specs) ]
        @test isapprox(mean(ss[1]), -0.05240393674216127, atol=1.0e-8)
        @test isapprox(mean(ss[2]), -0.02711953072380802, atol=1.0e-8)
        @test isapprox(mean(ss[3]), 0.00900794633309199, atol=1.0e-8)
        @test isapprox(mean(ss[4]), 0.1068253967760764, atol=1.0e-8)
        @test isapprox(mean(ss[5]), 0.11115156896689808, atol=1.0e-8)
        @test isapprox(std(ss[1]), 1.00731886928230457, atol=1.0e-8)
        @test isapprox(std(ss[2]), 1.0422450376085455, atol=1.0e-8)
        @test isapprox(std(ss[3]), 1.0294978144777747, atol=1.0e-8)
        @test isapprox(std(ss[4]), 0.9931443781679111, atol=1.0e-8)
        @test isapprox(std(ss[5]), 0.9858796014395503, atol=1.0e-8)

        n=fit(Normal, ss[1])
        @test isapprox(n.μ, -0.05240393674216127, atol = 1.0e-8)
        @test isapprox(n.σ, 1.007067008078397, atol = 1.0e-8)
    end

    @testset "uv" begin
        spec = loadspectrum(joinpath(@__DIR__,"K412 spectra","III-E K412[0][4].msa"))
        u = uv(spec)
        @test eltype(u)==UncertainValue
        @test length(u) == length(spec)
        @test all(value(u[i]) == spec[i] for i in eachindex(spec))
        @test all(value(u[i])==0.0 || σ(u[i]) == sqrt(spec[i]) for i in eachindex(spec))
        @test all(value(u[i])!=0.0 || σ(u[i]) == 1.0 for i in eachindex(spec))
    end
end
