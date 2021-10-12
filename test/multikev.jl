using NeXLSpectrum
using NeXLMatrixCorrection
using Test

@testset "15 keV vs 20 keV" begin    
    al2o3 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE Al2O3[$i][4].msa") for i in 0:4)))
    sio2 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE SiO2[$i][4].msa") for i in 0:4)))
    fe = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE Fe[$i][4].msa") for i in 0:4)))
    caf2 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE CaF2[$i][4].msa") for i in 0:4)))
    mgo = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE MgO[$i][4].msa") for i in 0:4)))
    unk20 = loadspectrum.(joinpath(@__DIR__,"Multi-keV", "20 keV", "IIIE K412[$i][4].msa") for i in 0:4)
    
    refs20 = references( [
        reference( n"Al", al2o3 ),
        reference( n"Si", sio2 ),
        reference( n"O", sio2 ),
        reference( n"Fe", fe ),
        reference( n"Ca", caf2 ),
        reference( n"Mg", mgo ),
    ], 132.0)

    q20 = map(s->quantify(s, refs20), unk20)
    # show(NeXLMatrixCorrection.describe(q20))
    m20 = mean(material.(q20))
    @test isapprox(value(m20[n"O"]), 0.4402, atol=0.0001)
    @test isapprox(value(m20[n"Mg"]), 0.1147, atol=0.0001)
    @test isapprox(value(m20[n"Al"]), 0.0484, atol=0.0001)
    @test isapprox(value(m20[n"Si"]), 0.2066, atol=0.0001)
    @test isapprox(value(m20[n"Ca"]), 0.1079, atol=0.0001)
    @test isapprox(value(m20[n"Fe"]), 0.0802, atol=0.0001)

    al2o3 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE Al2O3[$i][4].msa") for i in 0:4)))
    sio2 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE SiO2[$i][4].msa") for i in 0:4)))
    fe = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE Fe[$i][4].msa") for i in 0:4)))
    caf2 = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE CaF2[$i][4].msa") for i in 0:4)))
    mgo = sum(findsimilar(loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE MgO[$i][4].msa") for i in 0:4)))
    unk15 = loadspectrum.(joinpath(@__DIR__,"Multi-keV", "15 keV", "IIIE K412[$i][4].msa") for i in 0:4)
    
    refs15 = references( [
        reference( n"Al", al2o3 ),
        reference( n"Si", sio2 ),
        reference( n"O", sio2 ),
        reference( n"Fe", fe ),
        reference( n"Ca", caf2 ),
        reference( n"Mg", mgo ),
    ], 132.0)

    q15 = map(s->quantify(s, refs15), unk15)
    # show(NeXLMatrixCorrection.describe(q15))
    m15 = mean(material.(q15))
    @test isapprox(value(m15[n"O"]), 0.4566, atol=0.0001)
    @test isapprox(value(m15[n"Mg"]), 0.1157, atol=0.0001)
    @test isapprox(value(m15[n"Al"]), 0.0485, atol=0.0001)
    @test isapprox(value(m15[n"Si"]), 0.2102, atol=0.0001)
    @test isapprox(value(m15[n"Ca"]), 0.1082, atol=0.0001)
    @test isapprox(value(m15[n"Fe"]), 0.0803, atol=0.0001)



    q15_20 = map(s->quantify(s, refs20), unk15)
    #show(NeXLMatrixCorrection.describe(q15_20))
    m15_20 = mean(material.(q15_20))
    @test isapprox(value(m15_20[n"O"]), 0.4708, atol=0.0001)
    @test isapprox(value(m15_20[n"Mg"]), 0.1176, atol=0.0001)
    @test isapprox(value(m15_20[n"Al"]), 0.0488, atol=0.0001)
    @test isapprox(value(m15_20[n"Si"]), 0.2073, atol=0.0001)
    @test isapprox(value(m15_20[n"Ca"]), 0.1075, atol=0.0001)
    @test isapprox(value(m15_20[n"Fe"]), 0.0825, atol=0.0001)

    q20_15 = map(s->quantify(s, refs15), unk20)
    # show(NeXLMatrixCorrection.describe(q20_15))
    m20_15 = mean(material.(q20_15))
    @test isapprox(value(m20_15[n"O"]), 0.4278, atol=0.0001)
    @test isapprox(value(m20_15[n"Mg"]), 0.1128, atol=0.0001)
    @test isapprox(value(m20_15[n"Al"]), 0.0481, atol=0.0001)
    @test isapprox(value(m20_15[n"Si"]), 0.2094, atol=0.0001)
    @test isapprox(value(m20_15[n"Ca"]), 0.1087, atol=0.0001)
    @test isapprox(value(m20_15[n"Fe"]), 0.0781, atol=0.0001)
end


