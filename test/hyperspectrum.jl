using Test
using Pkg
using Pkg.Artifacts
using NeXLSpectrum
using DataFrames
using Downloads: download
using Statistics

# This is the path to the Artifacts.toml we will manipulate
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")

# Query the `Artifacts.toml` file for the hash bound to the name "rplraw"
# (returns `nothing` if no such binding exists)
rplraw_hash = artifact_hash("rplraw", artifacts_toml)

# If the name was not bound, or the hash it was bound to does not exist, create it!
if rplraw_hash == nothing || !artifact_exists(rplraw_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    rplraw_hash = create_artifact() do artifact_dir
        println("Artifact dir: $artifact_dir")
        url = "https://drive.google.com/uc?export=download&id=1C93kn9-EIXXMDPcqJ9E4Xt4j9qfs5eeX"
        tarball = if VERSION >= v"1.6.0-beta1.0"
            Downloads.download(url)
        else
            Base.download(url)
        end
        Pkg.probe_platform_engines!()
        Pkg.unpack(tarball, artifact_dir, verbose = true)
        rm(tarball)
    end
    # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
    # just overwrite with the new content-hash.  Unless the source files change, we do not expect
    # the content hash to change, so this should not cause unnecessary version control churn.
    bind_artifact!(artifacts_toml, "rplraw", rplraw_hash)
end

# Get the path of the rplraw dataset, either newly created or previously generated.
# this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
rrpath = artifact_path(rplraw_hash)

@testset "HyperSpectrum" begin
    raw = readrplraw(joinpath(rrpath, "map[15]"))

    @test size(raw) == (2048, 128, 128)
    @test eltype(raw) == UInt16
    @test raw[640, 64, 63] == 0x0006
    @test raw[25,21,45] == 0x000A
    @test raw[29, 121, 115] == 0x000F
    @test raw[27, 123, 45] == 0x000C
    
    les = LinearEnergyScale(0.0, 10.0)
    props = Dict{Symbol,Any}(
        :LiveTime => 0.004,
        :BeamEnergy => 20.0e3,
        :ProbeCurrent => 3.04,
        :Name => "Map[15]",
    )
    hs = HyperSpectrum(les, props, raw, fov = (0.128, 0.128))

    @test name(hs)=="Map[15]"
    @test size(hs) == (128, 128)
    @test eltype(hs) == Spectrum{UInt16}
    @test hs[64, 64] isa Spectrum
    @test hs[64, 64][640] == 0x0007
    @test ndims(hs) == 2

    @test channel(999.9, hs) == 100
    @test channel(1000.0, hs) == 100 + 1
    @test channel(1009.9, hs) == 100 + 1
    @test energy(100, hs) == 99.0 * 10.0
    @test rangeofenergies(100, hs) == (990.0, 1000.0)
    @test NeXLSpectrum.properties(hs)[:BeamEnergy] == hs[:BeamEnergy]

    mp = maxpixel(hs)
    @test mp isa Spectrum{UInt16}
    @test mp[640] == 0x0013
    @test indexofmaxpixel(hs, 640) == CartesianIndex(51, 64)
    @test hs[indexofmaxpixel(hs, 640)][640] == 0x0013
    @test hs[:BeamEnergy] == 20.0e3
    @test livetime(hs,(1,1)) == 0.004
    # @test all(sum(hs, map(ci->counts(hs,ci,640) .> 3, CartesianIndices(hs))) .<= sum(hs))

    @test length(hs) == 128 * 128
    @test ndims(hs) == 2
    @test size(hs) == (128, 128)
    @test size(hs, 1) == 128
    @test size(hs, 2) == 128
    @test axes(hs) == axes(hs.counts)[2:end]
    @test axes(hs, 1) == axes(hs.counts, 2)
    @test eachindex(hs) == eachindex(hs.counts[1, :, :])
    @test eachindex(hs) == Base.OneTo(128 * 128)
    @test stride(hs, 2) == 2048 * 128
    @test stride(hs, 1) == 2048
    @test strides(hs) == (2048, 2048 * 128)
    @test depth(hs) == 2048
    @test dose(hs) == hs[:ProbeCurrent] * mean(hs.livetime)
    @test isapprox(dose(hs), 3.04 * 0.004, atol=1.0e-6)

    @test counts(hs[22, 34]) == hs.counts[:, 22, 34]
    ci = CartesianIndex(22, 34)
    @test counts(hs[200], UInt16) == hs.counts[:, 72, 2]
    @test counts(hs[ci], UInt16) == hs.counts[:, 22, 34]
    specs = hs[[200, 201, 202]]
    specs2 = hs[200:202]
    @test length(specs) == 3
    @test length(specs2) == 3
    @test counts(specs[2]) == counts(specs[2])
    @test eltype(specs) == Spectrum{UInt16}
    @test specs[1][:Cartesian] == CartesianIndex(72, 2)
    @test specs[1] == hs[200]
    @test counts(hs[CartesianIndex(2, 72)]) == counts(hs[2, 72])
    shs = NeXLSpectrum.region(hs, 11:20, 45:2:64)
    @test size(shs, 1) == length(11:20)
    @test size(shs, 2) == length(45:2:64)
    @test counts(shs[1, 1]) == counts(hs[11, 45])
    @test counts(shs[2, 2]) == counts(hs[12, 47])
    hs1000 = hs[1:1000]
    @test length(hs1000) == 1000
    @test hs1000[1000] == hs[(1000-1)%128+1, (1000-1)÷128+1]
    @test hs[33, 78][110] == counts(hs, CartesianIndex(33, 78), 110)
    @test hs[64, 63][640] == 0x0006
    @test hs[21, 45][25] == 0x000A
    @test hs[121, 115][29] == 0x000F
    @test hs[123, 45][27] == 0x000C

    @test axisname(hs,1) == :Y
    @test axisname(hs,2) == :X
    @test axisvalue(hs, 1, 1) == -0.064
    @test axisvalue(hs, 1, 128) == 0.064
    @test isapprox(axisvalue(hs, 1, 28), -0.064 + (28-1)*0.128/(128-1), atol=1.0e-8)
    @test axisvalue(hs, 2, 1) == -0.064
    @test axisvalue(hs, 2, 128) == 0.064
    @test isapprox(axisvalue(hs, 1, 73), -0.064 + (73-1)*0.128/(128-1), atol=1.0e-8)

    @test first(axisrange(hs, 1))==-0.064
    @test last(axisrange(hs, 1))==0.064
    @test first(axisrange(hs, 2))==-0.064
    @test last(axisrange(hs, 2))==0.064
     
    raw = nothing
    GC.gc()
end

@testset "QQHyperspec" begin
    hs = HyperSpectrum(
        LinearEnergyScale(0.0, 10.0),
        Dict{Symbol,Any}(
            :LiveTime => 0.01,
            :BeamEnergy => 20.0e3,
            :ProbeCurrent => 1.0,
            :Name => "Map[15]",
        ),
        readrplraw(joinpath(rrpath, "map[15]")),
        fov = (0.128, 0.128)
    )
    ffp = references(
        [
            reference(n"C", joinpath(rrpath, "standards", "C std.msa"), mat"C"),
            reference(n"Fe", joinpath(rrpath, "standards", "Fe std.msa"), mat"Fe"),
            reference(n"S", joinpath(rrpath, "standards", "FeS2 std.msa"), mat"FeS2"),
            reference(n"O", joinpath(rrpath, "standards", "MgO std.msa"), mat"MgO"),
            reference(n"Si", joinpath(rrpath, "standards", "Si std.msa"), mat"Si"),
        ],
        132.0,
    )
    res = fit_spectrum(hs, ffp, mode = :Fast)
    @test res isa Vector{KRatios}
    # Test if :Fast quant equivalent to :Full quant for a couple of pixels
    res12_23 = fit_spectrum(hs[12, 23], ffp)

    matches(kr1::KRatio, kr2::KRatio) =
        kr1.element == kr2.element && brightest(kr1.lines) == brightest(kr2.lines)

    for krs in res
        kr1 = krs[12, 23]
        for kr2 in kratios(res12_23)
            if matches(kr1, kr2)
                @test equivalent(kr1.kratio, kr2.kratio)
            end
        end
    end

    res60_29 = fit_spectrum(hs[60, 29], ffp)
    for krs in res
        kr1 = krs[60, 29]
        for kr2 in kratios(res60_29)
            if matches(kr1, kr2)
                @test equivalent(kr1.kratio, kr2.kratio)
            end
        end
    end

    raw = nothing
    GC.gc()
end
