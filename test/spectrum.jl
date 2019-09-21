using Test
using NeXL

les = LinearEnergyScale(-495.0,5.0)

@testset "Energy Scale" begin
    @test energy(200,les)==-495.0+5.0*(200-1)
    @test energy(1,les)==-495.0+5.0*(1-1)
    @test energy(500,les)==-495.0+5.0*(500-1)

    @test channel(energy(200,les)+1.0,les)==200
    @test channel(energy(200,les)-1.0,les)==199
    @test channel(energy(200,les),les)==200

    @test channel(600.0, les)==220
    @test channel(1600.0, les)==round((1600.0+495.0)/5.0)+1
end

mnk = MnKaResolution(130.0)

@testset "Resolution" begin
    @test resolution(5898.7, mnk)==130.0
    @test isapprox(resolution(500.0, mnk),60.0,atol=1.0)
end

@testset "Spectrum" begin
    s1 = Spectrum(LinearEnergyScale(0.0,10.0),1:2048,
        Dict([(:BeamEnergy, 20.0), (:LiveTime,60.0)]))

    @test s1[:LiveTime]==60.0
    @test s1[100]==100
    @test s1[1000]==1000

    @test s1[100:10:150] ==  [ 100,110,120,130,140,150 ]
    s1[100] = 200
    s1[101:103] = [300, 400, 500]

    @test s1[99]==99
    @test s1[100]==200
    @test s1[101]==300
    @test s1[102]==400
    @test s1[103]==500
    @test s1[104]==104

    s1[:RealTime] = 90.0
    @test s1[:RealTime] == 90.0
end

@testset "Detector" begin
    eds = basicEDS(4096,5.0,-495.0,125.0)
    @test channelcount(eds) == 4096
    @test isapprox(resolution(10000.0,eds),160.0,atol=1.0)
    @test energy(2000,eds) == energy(2000,les)
    @test channel(energy(2000,eds),eds)==2000
end

@testset "readEMSA" begin
    sp = readEMSA("supplemental/K412 spectra/Al2O3_std.msa")
    @test sp[:BeamEnergy] == 20.0
    @test sp[:ProbeCurrent] == 1.10989
    @test sp[:Elevation] == 35.0
    @test sp[:LiveTime] == 1172.19288
    @test sp[:RealTime] == 1491.4828
    @test sp[:Owner] == "Unknown"
    @test sp[:Name] == "Al2O3 std"
    comp = sp[:Composition]

    @test comp[n"O"] == 47.0749
    @test comp[n"Al"] == 52.9251
    @test comp[n"Zr"] == 0.0
    @test density(comp) == 4.00

    chs = channel(1335.0,sp):channel(1688.0,sp)
    @test integrate(sp,chs) == 33143490

    low, high = 120:133, 172:180
    @test isapprox(integrate(sp,low,chs,high), 32097622,rtol=3.0e-3)

    ss = [ subsample(sp,0.1) for _ in 1:10 ]

    @test isapprox(sp[:LiveTime], sum(ss[i][:LiveTime] for i in eachindex(ss)), rtol=1.0e-9)
    @test isapprox(sp[:RealTime], sum(ss[i][:RealTime] for i in eachindex(ss)),rtol=1.0e-9)

    ch=channel(1450.0,sp)
    @test isapprox(sum(ss[i][ch] for i in eachindex(ss)),sp[ch],atol=4.0*sqrt(sp[ch])) # can fail occasionally...

end
