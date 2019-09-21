using Test
using NeXL

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

al2o3 = readEMSA("supplemental/K412 spectra/Al2O3_std.msa")
f2 = buildfilter(basicEDS(al2o3, 132.0))

alk = filter("Al K",al2o3,135:170,f2)
