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
        reference(n"Ti", joinpath(path,"Ti std.msa"), mat"Ti") ], det; filter=VariableWidthFilter)
    fs=fit_spectrum(k2496, refs)
    refs2 = references( [
        reference(n"Si", joinpath(path,"Si std.msa"), mat"Si"),
        reference(n"Ba", joinpath(path,"BaCl2 std.msa"), mat"BaCl2"),
        reference(n"O", joinpath(path,"MgO std.msa"), mat"MgO") ], det; filter=VariableWidthFilter);
    sanbornite = loadspectrum(joinpath(path,"Sanbornite std.msa"))
    sanbornite_std = fit_spectrum(sanbornite, refs2)
    fs_stds = standardize(fs, sanbornite_std, mat"BaSi2O5")
    kl = labels(fs_stds)
    ok = findfirst(l->n"O K-L3" in l.xrays, kl)
    @test isapprox( value(fs_stds.kratios, kl[ok]), 0.979407, atol=1.0e-5) 
    @test isapprox( σ(fs_stds.kratios, kl[ok]), 0.000818949, atol=1.0e-5) 
    bal = findfirst(l->n"Ba L3-M5" in l.xrays, kl)
    @test isapprox( value(fs_stds.kratios, kl[bal]), 0.821541, atol=1.0e-5) 
    @test isapprox( σ(fs_stds.kratios, kl[bal]), 0.00142848, atol=1.0e-5) 
    tik = findfirst(l->n"Ti K-L3" in l.xrays, kl)
    @test isapprox( value(fs_stds.kratios, kl[tik]), 0.022268, atol=1.0e-5) 
    @test isapprox( σ(fs_stds.kratios, kl[tik]), 0.000265, atol=1.0e-5) 
    q=material(quantify(fs_stds))
    @test isapprox(value(q[n"O"]), 0.3082074, atol=1.0e-5)
    @test isapprox(value(q[n"Si"]), 0.2219337, atol=1.0e-5)
    @test isapprox(value(q[n"Ti"]), 0.023002, atol=1.0e-5)
    @test isapprox(value(q[n"Ba"]), 0.42395893, atol=1.0e-5)

    stdks = mapreduce(append!,[n"Si",n"Ba",n"O"]) do elm
        extractStandards(sanbornite_std, elm, mat"BaSi2O5")
    end
    @assert all(isstandard, stdks)
    si_stds = extractStandards(sanbornite_std, n"Si", mat"BaSi2O5")

    @test all(isstandard, si_stds)
    @test length(si_stds)==1
    si_std = si_stds[1]
    @test element(si_std)==n"Si"
    @test isequal(si_std.standard, mat"Si")
    @test isequal(si_std.unkProps[:Composition], mat"BaSi2O5")
    @test isequal(si_std.stdProps[:Composition], mat"Si")

    # Standardize using k-ratios rather than a FitResult
    fs_stds2 = standardize(fs, stdks)
    @test isapprox( value(fs_stds2.kratios, kl[ok]), 0.979407, atol=1.0e-5) 
    @test isapprox( σ(fs_stds2.kratios, kl[ok]), 0.000818949, atol=1.0e-5) 
    @test isapprox( value(fs_stds2.kratios, kl[bal]), 0.8215410, atol=1.0e-5) 
    @test isapprox( σ(fs_stds2.kratios, kl[bal]), 0.00142848, atol=1.0e-5) 
    @test isapprox( value(fs_stds2.kratios, kl[tik]), 0.022268, atol=1.0e-5) 
    @test isapprox( σ(fs_stds2.kratios, kl[tik]), 0.000265, atol=1.0e-5) 

    ba_stds = extractStandards(sanbornite_std, [ n"Ba M5-N3", n"Ba L3-M5"], mat"BaSi2O5")
    @test all(isstandard, ba_stds)
    @test length(ba_stds)==2
    @test all(std->element(std)==n"Ba", ba_stds)
    @test all(std->isequal(std.standard,mat"BaCl2"), ba_stds)
    @test all(std->isequal(std.unkProps[:Composition],mat"BaSi2O5"), ba_stds)
end