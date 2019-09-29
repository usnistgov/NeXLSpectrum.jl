using Test
using NeXLCore
using NeXLSpectrum
using Printf

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

@testset "LLSQ" begin
    dir = @__DIR__
    path = "$(dir)/K412 spectra/"
    unks = readEMSA.(@sprintf("%sIII-E K412[%d][4].msa", path, i) for i in 0:4)
    al2o3 = readEMSA(path*"Al2O3 std.msa")
    caf2 = readEMSA(path*"CaF2 std.msa")
    fe = readEMSA(path*"Fe std.msa")
    mgo = readEMSA(path*"MgO std.msa")
    sio2 = readEMSA(path*"SiO2 std.msa");

    det = basicEDS(4096, 10.0, 0.0, 132.0)
    filt = buildfilter(det)

    ok = filter(sio2,34:66,ff)
    mgk = filter(mgo, 110:142, ff)
    alk = filter(al2o3,135:170,filt)
    sik = filter(sio2,159:196,ff)
    cak = filter(caf2, 345:422, ff)
    fel = filter(fe, 51:87, ff)
    feka = filter(fe , 615:666, ff)
    fekb = filter(fe, 677:735, ff)

    unk = filter(unks[1], ff)


end
