using NeXLSpectrum
using Test

@testset "K2496" begin
    path = joinpath(@__DIR__,"K2496")
    k2496 = loadspectrum(joinpath(path,"K2496_1.msa"))
    det = BasicEDS(4096, -480.40409, 5.00525, 132.0, 110, 
         Dict(KShell=>n"B", LShell=>n"Ca", MShell=>n"Cs"))
    apply(k2496,det)
    refs = references([
        reference(n"Si", joinpath(path,"Si std.msa"), mat"Si"),
        reference(n"Ba", joinpath(path,"BaCl2 std.msa"), mat"BaCl2"),
        reference(n"O", joinpath(path,"MgO std.msa"), mat"MgO"),
        reference(n"Ti", joinpath(path,"Ti std.msa"), mat"Ti") ], det)
    fs=fit_spectrum(k2496, refs)
    refs2 = references( [
        reference(n"Si", joinpath(path,"Si std.msa"), mat"Si"),
        reference(n"Ba", joinpath(path,"BaCl2 std.msa"), mat"BaCl2"),
        reference(n"O", joinpath(path,"MgO std.msa"), mat"MgO") ], det);
    sanbornite = loadspectrum(joinpath(path,"Sanbornite std.msa"))
    sanbornite_std = fit_spectrum(sanbornite, refs2)
    fs_stds = standardize(fs, sanbornite_std, mat"BaSi2O5")
    kl = labels(fs_stds)
    ok = findfirst(l->n"O K-L3" in l.xrays, kl)
    @test isapprox( value(kl[ok], fs_stds.kratios), 0.977767, atol=1.0e-6) 
    @test isapprox( σ(kl[ok], fs_stds.kratios), 0.000818949, atol=1.0e-6) 
    bal = findfirst(l->n"Ba L3-M5" in l.xrays, kl)
    @test isapprox( value(kl[bal], fs_stds.kratios), 0.846352, atol=1.0e-6) 
    @test isapprox( σ(kl[bal], fs_stds.kratios), 0.00142848, atol=1.0e-6) 
    tik = findfirst(l->n"Ti K-L3" in l.xrays, kl)
    @test isapprox( value(kl[tik], fs_stds.kratios), 0.0222705, atol=1.0e-6) 
    @test isapprox( σ(kl[tik], fs_stds.kratios), 0.000265, atol=1.0e-6) 
    q=quantify(fs_stds)
    @test isapprox(value(material(q)[n"O"]), 0.307589, atol=1.0e-6)
    @test isapprox(value(material(q)[n"Si"]), 0.221913, atol=1.0e-6)
    @test isapprox(value(material(q)[n"Ti"]), 0.022855, atol=1.0e-6)
    @test isapprox(value(material(q)[n"Ba"]), 0.423982, atol=1.0e-6)
end
