using NeXLSpectrum
using DataFrames
using Test


@testset "MultiDet" begin
    al2o3 = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE Al2O3[0]")
    sio2 = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE SiO2[0]")
    fe = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE Fe[0]")
    mgo = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE MgO[0]")
    caf2 = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE CaF2[0]")
    k412 = loadmultispec(joinpath(@__DIR__,"MultiSpec"),"IIIE K412[0]")


    comparesum(specs, sumspec) = isequal(multisum(specs).counts, sumspec.counts)
    @test comparesum(al2o3, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Al2O3[0][4].msa")))
    @test comparesum(sio2, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE SiO2[0][4].msa")))
    @test comparesum(fe, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Fe[0][4].msa")))
    @test comparesum(mgo, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE MgO[0][4].msa")))
    @test comparesum(caf2, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE CaF2[0][4].msa")))
    @test comparesum(k412, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE K412[0][4].msa")))

    @test isapprox(multisum(al2o3)[:LiveTime], loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Al2O3[0][4].msa"))[:LiveTime], atol=0.001)
    @test isapprox(multisum(sio2)[:LiveTime], loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE SiO2[0][4].msa"))[:LiveTime], atol=0.001)
    @test isapprox(multisum(fe)[:LiveTime], loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Fe[0][4].msa"))[:LiveTime], atol=0.001)
    @test isapprox(multisum(k412)[:LiveTime], loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE K412[0][4].msa"))[:LiveTime], atol=0.001)

    @test isapprox(dose(multisum(al2o3)), dose(loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Al2O3[0][4].msa"))), atol=0.001)
    @test isapprox(dose(multisum(sio2)), dose(loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE SiO2[0][4].msa"))), atol=0.001)
    @test isapprox(dose(multisum(fe)), dose(loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Fe[0][4].msa"))), atol=0.001)
    @test isapprox(dose(multisum(k412)), dose(loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE K412[0][4].msa"))), atol=0.001)

    comparemean(specs, sumspec) = isequal(multimean(specs).counts, sumspec.counts ./ length(specs))
    @test comparemean(al2o3, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Al2O3[0][4].msa")))
    @test comparemean(sio2, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE SiO2[0][4].msa")))
    @test comparemean(fe, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE Fe[0][4].msa")))
    @test comparemean(mgo, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE MgO[0][4].msa")))
    @test comparemean(caf2, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE CaF2[0][4].msa")))
    @test comparemean(k412, loadspectrum(joinpath(@__DIR__,"MultiSpec","IIIE K412[0][4].msa")))

    refs = references([
        reference(n"Si", multisum(sio2), mat"SiO2"),
        reference(n"Al", multisum(al2o3), mat"Al2O3"),
        reference(n"O", multisum(al2o3), mat"Al2O3"),
        reference(n"Fe", multisum(fe), mat"Fe"),
        reference(n"Mg", multisum(mgo), mat"MgO"),
        reference(n"Ca", multisum(caf2), mat"CaF2"),
    ], 132.0; filter=VariableWidthFilter )

    q = quantify(multisum(k412), refs)
    @test isapprox(value(q.comp[n"Al"]), 0.049, atol = 0.001)
    @test isapprox(value(q.comp[n"Ca"]), 0.108, atol = 0.001)
    @test isapprox(value(q.comp[n"Fe"]), 0.080, atol = 0.001)
    @test isapprox(value(q.comp[n"Mg"]), 0.117, atol = 0.001)
    @test isapprox(value(q.comp[n"Si"]), 0.206, atol = 0.001)
    @test isapprox(value(q.comp[n"O"]),  0.4571, atol = 0.001)

    @test all(isapprox(ms,v, atol=0.0001) for (ms, v) in zip(multiscore(al2o3), [0.0321, -0.0689, -0.0149, 0.0440 ]))
    @test all(isapprox(ms,v, atol=0.0001) for (ms, v) in zip(multiscore(sio2), [ 0.0478, -0.0473, -0.0287, 0.0196 ]))

    @test isapprox(multirank(k412), 0.0898, atol=0.0001)
    @test isapprox(multirank(fe), 0.0396, atol=0.0001)

end