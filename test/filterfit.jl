using Test
using NeXLCore
using NeXLSpectrum
using NeXLUncertainties
using Printf
using BenchmarkTools

@testset "Filter" begin
    eds = basicEDS(2048,10.0,0.0,135.0)
    filt = buildfilter(eds)
    # Each row sums to zero
    @test all(isapprox(sum(row),0.0,atol=1.0e-8) for row in eachrow(filt))
    # Symmetric about the center line
    @test all(isapprox(sum(filt[r,1:r-1]),sum(filt[r,r+1:end]),atol=1.0e-8) for r in 2:(size(filt)[1]-1))
    # Positive in the center
    @test all(row[r] ≥ 0.0 for (r,row) in enumerate(eachrow(filt)))
    # Symmetric one row off
    @test all( filt[r,r-1]==filt[r,r+1] for r in 2:(size(filt)[1]-1))
end

@testset "LLSQ_K412_1" begin
    dir = @__DIR__
    path = "$(dir)/K412 spectra/"
    unks = readEMSA.(@sprintf("%sIII-E K412[%d][4].msa", path, i) for i in 0:4)
    al2o3 = readEMSA(path*"Al2O3 std.msa")
    caf2 = readEMSA(path*"CaF2 std.msa")
    fe = readEMSA(path*"Fe std.msa")
    mgo = readEMSA(path*"MgO std.msa")
    sio2 = readEMSA(path*"SiO2 std.msa");

    det = basicEDS(4096, 10.0, 0.0, 132.0)
    ff = buildfilter(det)

    ok = filter(sio2,34:66,ff, 1.0/dose(sio2))
    mgk = filter(mgo, 110:142, ff, 1.0/dose(mgo))
    alk = filter(al2o3,135:170,ff, 1.0/dose(al2o3))
    sik = filter(sio2,159:196,ff, 1.0/dose(sio2))
    cak = filter(caf2, 345:422, ff, 1.0/dose(caf2))
    fel = filter(fe, 51:87, ff, 1.0/dose(fe))
    feka = filter(fe, 615:666, ff, 1.0/dose(fe))
    fekb = filter(fe, 677:735, ff, 1.0/dose(fe))

    fds = [ ok, mgk, alk, sik, cak, fel, feka, fekb ]

    unk = filter(unks[1], ff, 1.0/dose(unks[1]))

    ff = filterfit(unk, fds, fitcontiguousp)
    println("Performing the full generalized fit takes:")
    #@btime filterfit(unk, fds, fitcontiguousp)
    println("Performing the weighted fit takes:")
    #@btime filterfit(unk, fds, fitcontiguousw)

    @test isapprox(value(ok.identifier, ff), 0.6623, atol=0.0001)
    @test isapprox(value(fekb.identifier, ff), 0.0686, atol=0.0001)
    @test isapprox(value(mgk.identifier, ff), 0.1540, atol=0.0001)
    @test isapprox(value(alk.identifier, ff), 0.0691, atol=0.0001)
    @test isapprox(value(sik.identifier, ff), 0.3616, atol=0.0001)
    @test isapprox(value(cak.identifier, ff), 0.1952, atol=0.0001)
    @test isapprox(value(fel.identifier, ff), 0.0469, atol=0.0001)
    @test isapprox(value(feka.identifier, ff), 0.0715, atol=0.0001)

    @test isapprox(uncertainty(ok.identifier, ff), 0.000107, atol=0.000001)
    @test isapprox(uncertainty(mgk.identifier, ff), 3.810e-5, atol=1.0e-7)
    @test isapprox(uncertainty(alk.identifier, ff), 2.60e-5, atol=1.0e-7)
    @test isapprox(uncertainty(sik.identifier, ff), 6.21e-5, atol=1.0e-7)
    @test isapprox(uncertainty(cak.identifier, ff), 5.42e-5, atol=1.0e-7)
    @test isapprox(uncertainty(fel.identifier, ff), 0.0001809, atol=0.0001)
    @test isapprox(uncertainty(feka.identifier, ff), 4.36e-5, atol=1.0e-7)
    @test isapprox(uncertainty(fekb.identifier, ff), 0.0002178, atol=0.0001)
end

@testset "LLSQ_K412_2" begin
    dir = @__DIR__
    path = "$(dir)/K412 spectra/"
    unks = readEMSA.(@sprintf("%sIII-E K412[%d][4].msa", path, i) for i in 0:4)
    al2o3 = readEMSA(path*"Al2O3 std.msa")
    caf2 = readEMSA(path*"CaF2 std.msa")
    fe = readEMSA(path*"Fe std.msa")
    mgo = readEMSA(path*"MgO std.msa")
    sio2 = readEMSA(path*"SiO2 std.msa");

    det = basicEDSwICC(4096, 10.0, 0.0, 132.0)
    ff = buildfilter(det)

    ampl = 0.00005
    okroi = extent(characteristic(n"O",ktransitions),det,ampl)
    mgkroi = extent(characteristic(n"Mg",ktransitions),det,ampl)
    alkroi = extent(characteristic(n"Al",ktransitions),det,ampl)
    sikroi = extent(characteristic(n"Si",ktransitions),det,ampl)
    cakroi = extent(characteristic(n"Ca",ktransitions),det,ampl)
    felroi = extent(characteristic(n"Fe",ltransitions),det,ampl)
    fekaroi = extent(characteristic(n"Fe",kalpha),det,ampl)
    fekbroi = extent(characteristic(n"Fe",kother),det,ampl)

    ok = filter(sio2,okroi,ff, 1.0/dose(sio2))
    mgk = filter(mgo, mgkroi, ff, 1.0/dose(mgo))
    alk = filter(al2o3,alkroi,ff, 1.0/dose(al2o3))
    sik = filter(sio2,sikroi,ff, 1.0/dose(sio2))
    cak = filter(caf2, cakroi, ff, 1.0/dose(caf2))
    fel = filter(fe, felroi, ff, 1.0/dose(fe))
    feka = filter(fe, fekaroi, ff, 1.0/dose(fe))
    fekb = filter(fe, fekbroi, ff, 1.0/dose(fe))

    fds = [ ok, mgk, alk, sik, cak, fel, feka, fekb ]

    unk = filter(unks[1], ff, 1.0/dose(unks[1]))

    ff = filterfit(unk, fds, fitcontiguousp)
    println("Performing the full generalized fit takes:")
    #@btime filterfit(unk, fds, fitcontiguousp)
    println("Performing the weighted fit takes:")
    #@btime filterfit(unk, fds, fitcontiguousw)

    @test isapprox(value(ok.identifier, ff), 0.6623, atol=0.0001)
    @test isapprox(value(fekb.identifier, ff), 0.0686, atol=0.0001)
    @test isapprox(value(mgk.identifier, ff), 0.1540, atol=0.0001)
    @test isapprox(value(alk.identifier, ff), 0.0691, atol=0.0001)
    @test isapprox(value(sik.identifier, ff), 0.3616, atol=0.0001)
    @test isapprox(value(cak.identifier, ff), 0.1952, atol=0.0001)
    @test isapprox(value(fel.identifier, ff), 0.0469, atol=0.0001)
    @test isapprox(value(feka.identifier, ff), 0.0715, atol=0.0001)

    @test isapprox(uncertainty(ok.identifier, ff), 0.000107, atol=0.000001)
    @test isapprox(uncertainty(mgk.identifier, ff), 3.810e-5, atol=1.0e-7)
    @test isapprox(uncertainty(alk.identifier, ff), 2.60e-5, atol=1.0e-7)
    @test isapprox(uncertainty(sik.identifier, ff), 6.21e-5, atol=1.0e-7)
    @test isapprox(uncertainty(cak.identifier, ff), 5.42e-5, atol=1.0e-7)
    @test isapprox(uncertainty(fel.identifier, ff), 0.0001809, atol=0.0001)
    @test isapprox(uncertainty(feka.identifier, ff), 4.36e-5, atol=1.0e-7)
    @test isapprox(uncertainty(fekb.identifier, ff), 0.0002178, atol=0.0001)
end
