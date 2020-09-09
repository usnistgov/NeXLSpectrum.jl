using Test
using Pkg.Artifacts
using NeXLSpectrum
using DataFrames

# Download the necessary data using the Artifact mechanism from Google Drive
rrpath = artifact"rplraw"

@testset "HyperSpectrum" begin
    raw = readrplraw(joinpath(rrpath, "map[15]"))

    @test size(raw) == (2048, 128, 128)
    @test eltype(raw) == UInt16
    @test raw[640, 64, 64] == 0x0007

    les = LinearEnergyScale(0.0, 10.0)
    props = Dict{Symbol,Any}(:LiveTime => 0.004, :BeamEnergy => 20.0e3, :ProbeCurrent => 3.04, :Name=>"Map[15]")
    hs = HyperSpectrum(les, props, raw)

    @test size(hs) == (128, 128)
    @test eltype(hs) == Spectrum{UInt16}
    @test hs[64, 64] isa Spectrum
    @test hs[64, 64][640] == 0x0007
    @test ndims(hs)==2

    @test channel(999.9, hs)==100
    @test channel(1000.0, hs)==100+1
    @test channel(1009.9, hs)==100+1
    @test energy(100, hs) == 99.0*10.0
    @test rangeofenergies(100, hs) == (990.0,1000.0)
    @test NeXLSpectrum.properties(hs)[:BeamEnergy]==hs[:BeamEnergy]

    mp = maxpixel(hs)
    @test mp isa Spectrum{UInt16}
    @test mp[640] == 0x0013
    @test indexofmaxpixel(hs, 640) == CartesianIndex(64, 51)
    @test hs[indexofmaxpixel(hs, 640)][640] == 0x0013
    @test hs[:BeamEnergy]==20.0e3
    @test hs[:LiveTime]==0.004
    @test all(sum(hs, (sig, i) -> counts(hs, i, 640) > 3) .<= sum(hs))

    @test length(hs) == 128*128
    @test ndims(hs) == 2
    @test size(hs) == (128, 128)
    @test size(hs, 1) == 128
    @test size(hs, 2) == 128
    @test axes(hs) == axes(hs.counts)[2:end]
    @test axes(hs, 1) == axes(hs.counts, 2)
    @test eachindex(hs) == eachindex(hs.counts[1,:,:])
    @test eachindex(hs) == Base.OneTo(128*128)
    @test stride(hs, 2) == 2048*128
    @test stride(hs, 1) == 2048
    @test strides(hs) == (2048, 2048*128)
    @test depth(hs) == 2048
    @test dose(hs) == 3.04*0.004

    @test counts(hs[22,34]) == hs.counts[:, 22, 34]
    ci=CartesianIndex(22,34)
    @test counts(hs[200], UInt16)== hs.counts[:, 72, 2]
    @test counts(hs[ci], UInt16) == hs.counts[:, 22, 34]
    specs = hs[[200, 201, 202]]
    specs2 = hs[200:202]
    @test length(specs)==3
    @test length(specs2)==3
    @test counts(specs[2])==counts(specs[2])
    @test eltype(specs)==Spectrum{UInt16}
    @test specs[1][:Cartesian]==CartesianIndex(72, 2)
    @test specs[1]==hs[200]
    @test counts(hs[CartesianIndex(2,72)])==counts(hs[2,72])
    shs = NeXLSpectrum.region(hs,11:20,45:2:64)
    @test size(shs,1)==length(11:20)
    @test size(shs,2)==length(45:2:64)
    @test counts(shs[1,1])==counts(hs[11,45])
    @test counts(shs[2,2])==counts(hs[12,47])
    hs1000=hs[1:1000]
    @test length(hs1000)==1000
    @test hs1000[1000] == hs[(1000-1) % 128 + 1, (1000-1) ÷ 128 + 1]
    @test hs[33, 78][110]==counts(hs, CartesianIndex(33,78), 100)
    raw=nothing
    GC.gc()
end

@testset "QQHyperspec" begin
    les = LinearEnergyScale(0.0, 10.0)
    props = Dict{Symbol,Any}(:LiveTime => 0.01, :BeamEnergy => 20.0e3, :ProbeCurrent=>1.0, :Name=>"Map[15]")
    raw = readrplraw(joinpath(rrpath, "map[15]"))
    hs = HyperSpectrum(les, props, raw)
    mp = maxpixel(hs)
    cstd, festd, fes2std, mgostd, sistd = map(n->loadspectrum(joinpath(rrpath, "standards", "$n std.msa")), ("C", "Fe", "FeS2", "MgO", "Si"))
    det = matching(festd, 132.0, 10)
    refs = (
        ( cstd, n"C", mat"C" ),
        ( festd, n"Fe", mat"Fe" ),
        ( fes2std, n"S", mat"FeS2" ),
        ( mgostd, n"O", mat"MgO" ),
        ( sistd, n"Si", mat"Si" )
    )
    filt = buildfilter(det)
    frs = filterreferences(filt, refs...)
    qq = VectorQuant(frs, filt)
    res=fit(qq, hs) # Array{KRatios}
    @test res[1] isa KRatios
    @test all(map(a->isapprox(a..., atol=0.00001), zip((r.kratios[12,23] for r in res), ( 1.10457, 0.006895, 0.0, 0.0, 0.0, 0.063778, 0.0 ))))
    @test all(map(a->isapprox(a..., atol=0.00001), zip((r.kratios[60,29] for r in res), ( 1.070000, 0.01006, 0.36027, 0.20048, 0.025442, 1.11138, 0.002127 ))))
    raw = nothing
    GC.gc()
end
