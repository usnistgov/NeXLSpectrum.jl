using Test

@testset "HyperSpectrum" begin
    les = LinearEnergyScale(0.0,10.0)
    raw = readrplraw("C:\\Users\\nritchie\\Desktop\\map[15]",les,Dict{Symbol,Any}(:LiveTime=>0.004, :BeamEnergy=>20.0e3))

    @test size(raw)==(2048,128,128)
    @test eltype(raw)==UInt16
    @test raw[640,64,64]==0x0007

    mpr = maxpixel(raw)
    @test mpr isa Vector
    @test mpr[640]==0x0013
    @test raw[640, indexofmaxpixel(raw,640)]==0x0013

    hs = ashyperspectrum(raw)
    @test size(hs)==(128,128)
    @test eltype(hs)==Spectrum
    @test hs[64,64] isa Spectrum
    @test hs[64,64][640]==0x0007

    mp = maxpixel(hs)
    @test mp isa Spectrum
    @test mp[640]==0x0013
    @test indexofmaxpixel(hs, 640)==CartesianIndex(64, 51)
    @test hs[indexofmaxpixel(hs, 640)][640]==0x0013

    @test sum(raw)[:]==sum(hs)[:]
    @test all(sum(hs,(sig,i)->sig[640,i]>3) .<= sum(hs))
    @test all(sum(raw,(sig,i)->sig[640,i]>3) .<= sum(raw))
end;
