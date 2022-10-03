using Test
using NeXLSpectrum
using DataFrames
using Statistics

@testset "DirectFit" begin
    path = joinpath(@__DIR__,"K412 spectra")
    unks = loadspectrum.(joinpath(path, "III-E K412[$i][4].msa") for i = 0:4)
    al2o3 = loadspectrum(joinpath(path, "Al2O3 std.msa"))
    caf2 = loadspectrum(joinpath(path, "CaF2 std.msa"))
    fe = loadspectrum(joinpath(path, "Fe std.msa"))
    mgo = loadspectrum(joinpath(path, "MgO std.msa"))
    sio2 = loadspectrum(joinpath(path, "SiO2 std.msa"));

    det = matching(unks[1], 132.0)
    resp = detectorresponse(det, SDDEfficiency(ModeledWindow(MoxtekAP33())));

    drefs = references( [
        direct(n"Al", al2o3, mat"Al2O3"),
        direct(n"Ca", caf2, mat"CaF2"),
        direct(n"Fe", fe, mat"Fe"),
        direct(n"Mg", mgo, mat"MgO"),
        direct(n"Si", sio2, mat"SiO2"),
        direct(n"O", sio2, mat"SiO2"),
    ] , det, resp)

    dffrs = map(sp->fit_spectrum(sp, drefs), unks)

    df = asa(DataFrame, dffrs)
    @test isapprox(mean(df[:, "k[O K-L3 + 1 other, SiO2]"]), 0.638220, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Fe L3-M5 + 13 others, Fe]"]), 0.027000, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Mg K-L3 + 1 other, MgO]"]), 0.147659, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Al K-L3 + 3 others, Al2O3]"]), 0.068047, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Si K-L3 + 3 others, SiO2]"]), 0.353672, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Ca K-L3 + 3 others, CaF2]"]), 0.191788, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Fe K-L3 + 1 other, Fe]"]), 0.065704, atol=1.0e-6)
    @test isapprox(mean(df[:, "k[Fe K-M3 + 3 others, Fe]"]), 0.064218, atol=1.0e-6)

    qs = quantify.(dffrs)
    qdf = asa(DataFrame, qs)
    @test isapprox(mean(qdf[:, "O"]), 0.434617, atol=1.0e-6)
    @test isapprox(mean(qdf[:, "Mg"]), 0.115140, atol=1.0e-6)
    @test isapprox(mean(qdf[:, "Al"]), 0.0489210, atol=1.0e-6)
    @test isapprox(mean(qdf[:, "Si"]), 0.209054, atol=1.0e-6)
    @test isapprox(mean(qdf[:, "Ca"]), 0.108645, atol=1.0e-6)
    @test isapprox(mean(qdf[:, "Fe"]), 0.079175, atol=1.0e-6)

    # Compare to nominal K412
    @test isapprox(mean(qdf[:, "O"]), 0.427576, atol=0.008)
    @test isapprox(mean(qdf[:, "Mg"]), 0.116568, atol=2.0e-3)
    @test isapprox(mean(qdf[:, "Al"]), 0.0490621, atol=1.0e-3)
    @test isapprox(mean(qdf[:, "Si"]), 0.211983, atol=3.0e-3)
    @test isapprox(mean(qdf[:, "Ca"]), 0.108991, atol=1.0e-3)
    @test isapprox(mean(qdf[:, "Fe"]), 0.0774201, atol=2.2e-3)
end