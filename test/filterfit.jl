using Test
using NeXLSpectrum
using Statistics
using Printf

@testset "Filter Fitting" begin
    @testset "Filter" begin
        eds = simpleEDS(2048, 10.0, 0.0, 135.0)
        filt = buildfilter(eds)
        # Each row sums to zero
        @test all(isapprox(sum(NeXLSpectrum.filterdata(filt, row)), 0.0, atol = 1.0e-8) for row in size(filt)[1])
        # Symmetric about the center line
        @test all(
            isapprox(
                sum(NeXLSpectrum.filterdata(filt, r)[1:r-1]),
                sum(NeXLSpectrum.filterdata(filt, r)[r+1:end]),
                atol = 1.0e-8,
            ) for r = 2:(size(filt)[1]-1)
        )
        # Positive in the center
        @test all(NeXLSpectrum.filterdata(filt, r)[r] ≥ 0.0 for r = 1:size(filt)[1])
        # Symmetric one row off
        @test all(
            NeXLSpectrum.filterdata(filt, r)[r-1] == NeXLSpectrum.filterdata(filt, r)[r+1] for r = 2:(size(filt)[1]-1)
        )
    end

    @testset "LLSQ_K412_1" begin
        dir = @__DIR__
        path = "$(dir)/K412 spectra/"
        unks = readEMSA.(@sprintf("%sIII-E K412[%d][4].msa", path, i) for i = 0:4)
        al2o3 = readEMSA(path * "Al2O3 std.msa")
        caf2 = readEMSA(path * "CaF2 std.msa")
        fe = readEMSA(path * "Fe std.msa")
        mgo = readEMSA(path * "MgO std.msa")
        sio2 = readEMSA(path * "SiO2 std.msa")

        det = simpleEDS(4096, 10.0, 0.0, 132.0)
        ff = buildfilter(det)

        ok = tophatfilter(CharXRayLabel(sio2, 34:66, characteristic(n"O", ktransitions)), ff, 1.0 / dose(sio2))
        mgk = tophatfilter(CharXRayLabel(mgo, 110:142, characteristic(n"Mg", ktransitions)), ff, 1.0 / dose(mgo))
        alk = tophatfilter(CharXRayLabel(al2o3, 135:170, characteristic(n"Al", ktransitions)), ff, 1.0 / dose(al2o3))
        sik = tophatfilter(CharXRayLabel(sio2, 159:196, characteristic(n"Si", ktransitions)), ff, 1.0 / dose(sio2))
        cak = tophatfilter(CharXRayLabel(caf2, 345:422, characteristic(n"Ca", ktransitions)), ff, 1.0 / dose(caf2))
        fel = tophatfilter(CharXRayLabel(fe, 51:87, characteristic(n"Fe", ltransitions)), ff, 1.0 / dose(fe))
        feka = tophatfilter(CharXRayLabel(fe, 615:666, characteristic(n"Fe", kalpha)), ff, 1.0 / dose(fe))
        fekb = tophatfilter(CharXRayLabel(fe, 677:735, characteristic(n"Fe", kbeta)), ff, 1.0 / dose(fe))

        fds = [ok, mgk, alk, sik, cak, fel, feka, fekb]

        unk = tophatfilter(unks[1], ff, 1.0 / dose(unks[1]))

        #println("Performing the weighted fit takes:")
        ff = filterfit(unk, fds)
        #@btime filterfit(unk, fds)
        #println("Performing the full generalized fit takes:")
        #@btime filterfit(unk, fds, fitcontiguousp)
        #@btime filterfit(unk, fds, fitcontiguousw)

        ## The comparison is against the k-ratios from DTSA-II.
        # The results won't be identical because the filters and other assumptions are different.
        @test isapprox(value(ok.identifier, ff), 0.6529, atol = 0.0016)
        @test isapprox(value(fekb.identifier, ff), 0.0665, atol = 0.0002)
        @test isapprox(value(mgk.identifier, ff), 0.1473, atol = 0.0004)
        @test isapprox(value(alk.identifier, ff), 0.0668, atol = 0.0005)
        @test isapprox(value(sik.identifier, ff), 0.3506, atol = 0.0008)
        @test isapprox(value(cak.identifier, ff), 0.1921, atol = 0.0001)
        @test isapprox(value(fel.identifier, ff), 0.0418, atol = 0.0002)
        @test isapprox(value(feka.identifier, ff), 0.0669, atol = 0.0001)

        @test isapprox(uncertainty(ok.identifier, ff), 0.00061, atol = 0.0001)
        @test isapprox(uncertainty(mgk.identifier, ff), 0.00013, atol = 0.00005)
        @test isapprox(uncertainty(alk.identifier, ff), 0.00012, atol = 0.00005)
        @test isapprox(uncertainty(sik.identifier, ff), 0.00024, atol = 0.00005)
        @test isapprox(uncertainty(cak.identifier, ff), 0.00022, atol = 0.00005)
        @test isapprox(uncertainty(fel.identifier, ff), 0.00031, atol = 0.00006)
        @test isapprox(uncertainty(feka.identifier, ff), 0.00019, atol = 0.00005)
        @test isapprox(uncertainty(fekb.identifier, ff), 0.00078, atol = 0.0002)
    end;

    @testset "LLSQ_K412_2" begin
        dir = @__DIR__
        path = "$(dir)/K412 spectra/"
        unks = readEMSA.(@sprintf("%sIII-E K412[%d][4].msa", path, i) for i = 0:4)
        al2o3 = readEMSA(path * "Al2O3 std.msa")
        caf2 = readEMSA(path * "CaF2 std.msa")
        fe = readEMSA(path * "Fe std.msa")
        mgo = readEMSA(path * "MgO std.msa")
        sio2 = readEMSA(path * "SiO2 std.msa")

        det = BasicEDS(4096, 0.0, 10.0, 132.0, 10, Dict(Shell(1) => n"Be", Shell(2) => n"Sc", Shell(3) => n"Cs", Shell(4) => n"Pu"))
        ff = buildfilter(det)
        e0 = sameproperty(unks, :BeamEnergy)

        ampl = 0.00005
        oroi = charXRayLabels(sio2, n"O", [n"Si", n"O"], det, e0, ampl=ampl)
        siroi = charXRayLabels(sio2, n"Si", [n"Si", n"O"], det, e0, ampl=ampl)
        mgroi = charXRayLabels(mgo, n"Mg", [n"Mg", n"O"], det, e0, ampl=ampl)
        alroi = charXRayLabels(al2o3, n"Al", [n"Al", n"O"], det, e0, ampl=ampl)
        caroi = charXRayLabels(caf2, n"Ca", [n"Ca", n"F"], det, e0, ampl=ampl)
        @test length(caroi) == 1
        feroi = charXRayLabels(fe, n"Fe", [n"Fe"], det, e0, ampl=ampl)
        @test length(feroi) == 3

        ok = tophatfilter(oroi, ff, 1.0 / dose(sio2))
        mgk = tophatfilter(mgroi, ff, 1.0 / dose(mgo))
        alk = tophatfilter(alroi, ff, 1.0 / dose(al2o3))
        sik = tophatfilter(siroi, ff, 1.0 / dose(sio2))
        cak = tophatfilter(caroi, ff, 1.0 / dose(caf2))
        fekl = tophatfilter(feroi, ff, 1.0 / dose(fe))

        fds = collect(Iterators.flatten((ok, mgk, alk, sik, cak, fekl)))

        unk = tophatfilter(unks[1], ff, 1.0 / dose(unks[1]))

        ff = filterfit(unk, fds)
        #println("Performing the full generalized fit takes:")
        #@btime filterfit(unk, fds)
        #println("Performing the weighted fit takes:")
        #@btime filterfit(unk, fds, fitcontiguousw)

        @test isapprox(value(oroi[1], ff), 0.6624, atol = 0.0006)
        @test isapprox(value(mgroi[1], ff), 0.14728, atol = 0.0007)
        @test isapprox(value(alroi[1], ff), 0.06679, atol = 0.0006)
        @test isapprox(value(siroi[1], ff), 0.35063, atol = 0.0009)
        @test isapprox(value(caroi[1], ff), 0.19213, atol = 0.0003)
        @test isapprox(value(feroi[1], ff), 0.04185, atol = 0.0004)
        @test isapprox(value(feroi[2], ff), 0.06693, atol = 0.0001)
        @test isapprox(value(feroi[3], ff), 0.06652, atol = 0.0007)

        @test isapprox(uncertainty(oroi[1], ff), 0.00061, atol = 0.0001)
        @test isapprox(uncertainty(mgroi[1], ff), 0.00013, atol = 0.00004)
        @test isapprox(uncertainty(alroi[1], ff), 0.00012, atol = 0.00003)
        @test isapprox(uncertainty(siroi[1], ff), 0.00024, atol = 0.00003)
        @test isapprox(uncertainty(caroi[1], ff), 0.00022, atol = 0.00001)
        @test isapprox(uncertainty(feroi[1], ff), 0.00031, atol = 0.0001)
        @test isapprox(uncertainty(feroi[2], ff), 0.00019, atol = 0.00005)
        @test isapprox(uncertainty(feroi[3], ff), 0.00078, atol = 0.0002)
    end

    @testset "ADM6005a" begin
        dir = @__DIR__
        path = "$(dir)/ADM6005a spectra/"
        unks = readEMSA.("$(path)ADM-6005a_$(i).msa" for i = 1:15)
        al = readEMSA("$(path)Al std.msa")
        caf2 = readEMSA("$(path)CaF2 std.msa")
        fe = readEMSA("$(path)Fe std.msa")
        ge = readEMSA("$(path)Ge std.msa")
        si = readEMSA("$(path)Si std.msa")
        sio2 = readEMSA("$(path)SiO2 std.msa")
        ti = readEMSA("$(path)Ti trimmed.msa")
        zn = readEMSA("$(path)Zn std.msa")

        det = matching(unks[1], 128.0, 110, Dict(Shell(1) => n"Be", Shell(2) => n"Sc", Shell(3) => n"Cs", Shell(4) => n"Pu"))
        ff = buildfilter(det)

        ampl = 1e-4
        e0 = sameproperty(unks, :BeamEnergy)
        alroi = charXRayLabels(al, n"Al", [n"Al"], det, e0, ampl=ampl)
        caroi = charXRayLabels(caf2, n"Ca", [n"Ca",n"F"], det, e0, ampl=ampl)
        feroi = charXRayLabels(fe, n"Fe", [n"Fe"], det, e0, ampl=ampl)
        geroi = charXRayLabels(ge, n"Ge", [n"Ge"], det, e0, ampl=ampl)
        oroi = charXRayLabels(sio2, n"O", [n"Si",n"O"], det, e0, ampl=ampl)
        siroi = charXRayLabels(si, n"Si", [n"Si"], det, e0, ampl=ampl)
        tiroi = charXRayLabels(ti, n"Ti", [n"Ti"], det, e0, ampl=ampl)
        znroi = charXRayLabels(zn, n"Zn", [n"Zn"], det, e0, ampl=ampl)

        alk = tophatfilter(alroi, ff, 1.0 / dose(al))
        cak = tophatfilter(caroi, ff, 1.0 / dose(caf2))
        felk = tophatfilter(feroi, ff, 1.0 / dose(fe))
        gelk = tophatfilter(geroi, ff, 1.0 / dose(ge))
        ok = tophatfilter(oroi, ff, 1.0 / dose(sio2))
        sik = tophatfilter(siroi, ff, 1.0 / dose(si))
        tilk = tophatfilter(tiroi, ff, 1.0 / dose(ti))
        znlk = tophatfilter(znroi, ff, 1.0 / dose(zn))

        fds = collect(Iterators.flatten( ( alk, cak, felk, gelk, ok, sik, tilk, znlk ) ))

        res = FilterFitResult[]
        for i = 1:15
            unk = tophatfilter(unks[i], ff, 1.0 / dose(unks[i]))
            push!(res, filterfit(unk, fds))
        end

        # Compare against DTSA-II values
        @test isapprox(mean(values(oroi[1], res)), 0.4923, rtol=0.003)
        @test isapprox(mean(values(siroi[1], res)), 0.0214, atol=0.013)
        @test isapprox(mean(values(alroi[1], res)), 0.0281, rtol=0.004)
        @test isapprox(mean(values(caroi[1], res)), 0.1211, rtol=0.0025)
        @test isapprox(mean(values(znroi[1], res)), 0.0700, rtol=0.05)
        @test isapprox(mean(values(znroi[2], res)), 0.1115, atol=0.0005)
        @test isapprox(mean(values(znroi[3], res)), 0.1231, rtol=0.01)
        @test isapprox(mean(values(tiroi[1], res)), 0.0541, rtol=0.26)
        @test isapprox(mean(values(tiroi[2], res)), 0.064, rtol=0.0002)
        @test isapprox(mean(values(tiroi[3], res)), 0.064, rtol=0.06)
        @test isapprox(mean(values(feroi[1], res)), 0.0, atol=0.001)
        @test isapprox(mean(values(feroi[2], res)), 0.0, atol=0.0004)
        @test isapprox(mean(values(feroi[3], res)), 0.0, atol=0.001)
        @test isapprox(mean(values(geroi[1], res)), 0.1789, rtol=0.01)
        @test isapprox(mean(values(geroi[2], res)), 0.2628, atol=0.001)
        @test isapprox(mean(values(geroi[3], res)), 0.279, atol=0.011)
    end
end
