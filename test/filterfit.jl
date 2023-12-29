using Test
using NeXLSpectrum
using Statistics
using LinearAlgebra
using DataFrames

@testset "Filter Fitting" begin
    @testset "Filter" begin
        eds = simpleEDS(2048, 10.0, 0.0, 135.0)
        filt = buildfilter(eds)
        # Each row sums to zero
        @test all(
            isapprox(sum(NeXLSpectrum.filterdata(filt, row)), 0.0, atol = 1.0e-8) for
            row in size(filt)[1]
        )
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
            NeXLSpectrum.filterdata(filt, r)[r-1] == NeXLSpectrum.filterdata(filt, r)[r+1]
            for r = 2:(size(filt)[1]-1)
        )
        # Check the old and new ways are equivalent
        @test NeXLSpectrum.filterdata(filt, 1:size(filt)[1]) == NeXLSpectrum.filterdata(filt)
    end

    @testset "LLSQ_K412_1" begin
        path = joinpath(@__DIR__, "K412 spectra")
        unks = loadspectrum.(joinpath(path, "III-E K412[$i][4].msa") for i = 0:4)
        al2o3 = loadspectrum(joinpath(path, "Al2O3 std.msa"))
        caf2 = loadspectrum(joinpath(path, "CaF2 std.msa"))
        fe = loadspectrum(joinpath(path, "Fe std.msa"))
        mgo = loadspectrum(joinpath(path, "MgO std.msa"))
        sio2 = loadspectrum(joinpath(path, "SiO2 std.msa"))

        det = simpleEDS(4096, 10.0, 0.0, 132.0)
        ff = buildfilter(det)

        ok = tophatfilter(
            CharXRayLabel(sio2, 34:66, characteristic(n"O", ktransitions)),
            ff,
            1.0 / dose(sio2),
        )
        mgk = tophatfilter(
            CharXRayLabel(mgo, 110:142, characteristic(n"Mg", ktransitions)),
            ff,
            1.0 / dose(mgo),
        )
        alk = tophatfilter(
            CharXRayLabel(al2o3, 135:170, characteristic(n"Al", ktransitions)),
            ff,
            1.0 / dose(al2o3),
        )
        sik = tophatfilter(
            CharXRayLabel(sio2, 159:196, characteristic(n"Si", ktransitions)),
            ff,
            1.0 / dose(sio2),
        )
        cak = tophatfilter(
            CharXRayLabel(caf2, 345:422, characteristic(n"Ca", ktransitions)),
            ff,
            1.0 / dose(caf2),
        )
        fel = tophatfilter(
            CharXRayLabel(fe, 51:87, characteristic(n"Fe", ltransitions)),
            ff,
            1.0 / dose(fe),
        )
        feka = tophatfilter(
            CharXRayLabel(fe, 615:666, characteristic(n"Fe", kalpha)),
            ff,
            1.0 / dose(fe),
        )
        fekb = tophatfilter(
            CharXRayLabel(fe, 677:735, characteristic(n"Fe", kbeta)),
            ff,
            1.0 / dose(fe),
        )

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
        @test isapprox(NeXLUncertainties.value(ff, ok.label), 0.6529, atol = 0.0016)
        @test isapprox(NeXLUncertainties.value(ff, fekb.label), 0.0665, atol = 0.0002)
        @test isapprox(NeXLUncertainties.value(ff, mgk.label), 0.1473, atol = 0.0004)
        @test isapprox(NeXLUncertainties.value(ff, alk.label), 0.0668, atol = 0.0005)
        @test isapprox(NeXLUncertainties.value(ff, sik.label), 0.3506, atol = 0.0008)
        @test isapprox(NeXLUncertainties.value(ff, cak.label), 0.1921, atol = 0.0001)
        @test isapprox(NeXLUncertainties.value(ff, fel.label), 0.0418, atol = 0.0002)
        @test isapprox(NeXLUncertainties.value(ff, feka.label), 0.0669, atol = 0.0001)

        @test isapprox(σ(ff, ok.label), 0.00081, atol = 0.0001)
        @test isapprox(σ(ff, mgk.label), 0.00018, atol = 0.00005)
        @test isapprox(σ(ff, alk.label), 0.00012, atol = 0.00005)
        @test isapprox(σ(ff, sik.label), 0.00024, atol = 0.00005)
        @test isapprox(σ(ff, cak.label), 0.00022, atol = 0.00005)
        @test isapprox(σ(ff, fel.label), 0.00043, atol = 0.00006)
        @test isapprox(σ(ff, feka.label), 0.00019, atol = 0.00005)
        @test isapprox(σ(ff, fekb.label), 0.00078, atol = 0.0002)
    end

    @testset "LLSQ_K412_2" begin
        path = joinpath(@__DIR__, "K412 spectra")
        unks = loadspectrum.(joinpath(path, "III-E K412[$i][4].msa") for i = 0:4)
        al2o3 = loadspectrum(joinpath(path, "Al2O3 std.msa"))
        caf2 = loadspectrum(joinpath(path, "CaF2 std.msa"))
        fe = loadspectrum(joinpath(path, "Fe std.msa"))
        mgo = loadspectrum(joinpath(path, "MgO std.msa"))
        sio2 = loadspectrum(joinpath(path, "SiO2 std.msa"))

        det = BasicEDS(
            4096,
            0.0,
            10.0,
            132.0,
            10,
            Dict(
                Shell(1) => n"Be",
                Shell(2) => n"Sc",
                Shell(3) => n"Cs",
                Shell(4) => n"Pu",
            ),
        )
        ff = buildfilter(det)
        e0 = sameproperty(unks, :BeamEnergy)

        ampl = 0.00005
        oroi = charXRayLabels(sio2, n"O", Set( ( n"Si", n"O" ) ), det, e0, ampl = ampl)
        siroi = charXRayLabels(sio2, n"Si", Set( (n"Si", n"O") ), det, e0, ampl = ampl)
        mgroi = charXRayLabels(mgo, n"Mg", Set( (n"Mg", n"O") ), det, e0, ampl = ampl)
        alroi = charXRayLabels(al2o3, n"Al", Set( (n"Al", n"O") ), det, e0, ampl = ampl)
        caroi = charXRayLabels(caf2, n"Ca", Set( (n"Ca", n"F") ), det, e0, ampl = ampl)
        @test length(caroi) == 1
        feroi = charXRayLabels(fe, n"Fe", Set( (n"Fe", ) ), det, e0, ampl = ampl)
        @test length(feroi) == 3

        ok = tophatfilter(oroi, ff, 1.0 / dose(sio2))
        mgk = tophatfilter(mgroi, ff, 1.0 / dose(mgo))
        alk = tophatfilter(alroi, ff, 1.0 / dose(al2o3))
        sik = tophatfilter(siroi, ff, 1.0 / dose(sio2))
        cak = tophatfilter(caroi, ff, 1.0 / dose(caf2))
        fekl = tophatfilter(feroi, ff, 1.0 / dose(fe))

        fds = reduce(append!, ( ok, mgk, alk, sik, cak, fekl ))

        unk = tophatfilter(unks[1], ff, 1.0 / dose(unks[1]))

        ff = filterfit(unk, fds)
        #println("Performing the full generalized fit takes:")
        #@btime filterfit(unk, fds)
        #println("Performing the weighted fit takes:")
        #@btime filterfit(unk, fds, fitcontiguousw)

        @test isapprox(NeXLUncertainties.value(ff, oroi[1]), 0.6624, atol = 0.0001)
        @test isapprox(NeXLUncertainties.value(ff, mgroi[1]), 0.14728, atol = 0.0007)
        @test isapprox(NeXLUncertainties.value(ff, alroi[1]), 0.06679, atol = 0.0006)
        @test isapprox(NeXLUncertainties.value(ff, siroi[1]), 0.35063, atol = 0.0009)
        @test isapprox(NeXLUncertainties.value(ff, caroi[1]), 0.19213, atol = 0.0003)
        @test isapprox(NeXLUncertainties.value(ff, feroi[1]), 0.04185, atol = 0.0004)
        @test isapprox(NeXLUncertainties.value(ff, feroi[2]), 0.06693, atol = 0.0001)
        @test isapprox(NeXLUncertainties.value(ff, feroi[3]), 0.06652, atol = 0.0007)

        @test isapprox(σ(ff, oroi[1]), 0.00082, atol = 0.0001)
        @test isapprox(σ(ff, mgroi[1]), 0.00018, atol = 0.00004)
        @test isapprox(σ(ff, alroi[1]), 0.00016, atol = 0.00003)
        @test isapprox(σ(ff, siroi[1]), 0.00029, atol = 0.00003)
        @test isapprox(σ(ff, caroi[1]), 0.00023, atol = 0.00001)
        @test isapprox(σ(ff, feroi[1]), 0.00044, atol = 0.0001)
        @test isapprox(σ(ff, feroi[2]), 0.00016, atol = 0.00001)
        @test isapprox(σ(ff, feroi[3]), 0.00078, atol = 0.0002)

        # Compare to naive peak integration
        fekkr = kratio(unks[1], fe, 593:613, 636:647, 669:690)
        @test isapprox(
            NeXLUncertainties.value(ff, feroi[2]),
            NeXLUncertainties.value(fekkr),
            atol = 0.0005,
        )
        @test isapprox(σ(ff, feroi[2]), σ(fekkr), atol = 0.00004)

        cakkr = kratio(unks[1], caf2, 334:347, 365:375, 422:439)
        @test isapprox(
            NeXLUncertainties.value(ff, caroi[1]),
            NeXLUncertainties.value(cakkr),
            atol = 0.0008,
        )
        @test isapprox(σ(ff, caroi[1]), σ(cakkr), atol = 0.00007)
    end

    @testset "ADM6005a" begin
        path = joinpath(@__DIR__, "ADM6005a spectra")
        unks = loadspectrum.(joinpath(path, "ADM-6005a_$(i).msa") for i = 1:15)
        al = loadspectrum(joinpath(path, "Al std.msa"))
        caf2 = loadspectrum(joinpath(path, "CaF2 std.msa"))
        fe = loadspectrum(joinpath(path, "Fe std.msa"))
        ge = loadspectrum(joinpath(path, "Ge std.msa"))
        si = loadspectrum(joinpath(path, "Si std.msa"))
        sio2 = loadspectrum(joinpath(path, "SiO2 std.msa"))
        ti = loadspectrum(joinpath(path, "Ti trimmed.msa"))
        zn = loadspectrum(joinpath(path, "Zn std.msa"))

        det = matching(
            unks[1],
            128.0,
            110,
            Dict(
                Shell(1) => n"Be",
                Shell(2) => n"Sc",
                Shell(3) => n"Cs",
                Shell(4) => n"Pu",
            ),
        )
        ff = buildfilter(det)

        ampl = 1e-4
        e0 = sameproperty(unks, :BeamEnergy)
        alroi = charXRayLabels(al, n"Al", Set( ( n"Al", )), det, e0, ampl = ampl)
        caroi = charXRayLabels(caf2, n"Ca", [ n"Ca", n"F" ], det, e0, ampl = ampl)
        feroi = charXRayLabels(fe, n"Fe", Set( ( n"Fe", )), det, e0, ampl = ampl)
        geroi = charXRayLabels(ge, n"Ge", Set( ( n"Ge", )), det, e0, ampl = ampl)
        oroi = charXRayLabels(sio2, n"O", Set( ( n"Si", n"O", )), det, e0, ampl = ampl)
        siroi = charXRayLabels(si, n"Si", Set( ( n"Si", )), det, e0, ampl = ampl)
        tiroi = charXRayLabels(ti, n"Ti", Set( ( n"Ti", )), det, e0, ampl = ampl)
        znroi = charXRayLabels(zn, n"Zn", [ n"Zn" ], det, e0, ampl = ampl)

        alk = tophatfilter(alroi, ff, 1.0 / dose(al))
        cak = tophatfilter(caroi, ff, 1.0 / dose(caf2))
        felk = tophatfilter(feroi, ff, 1.0 / dose(fe))
        gelk = tophatfilter(geroi, ff, 1.0 / dose(ge))
        ok = tophatfilter(oroi, ff, 1.0 / dose(sio2))
        sik = tophatfilter(siroi, ff, 1.0 / dose(si))
        tilk = tophatfilter(tiroi, ff, 1.0 / dose(ti))
        znlk = tophatfilter(znroi, ff, 1.0 / dose(zn))

        fds = reduce(append!, ( alk, cak, felk, gelk, ok, sik, tilk, znlk ))

        res = FilterFitResult[]
        for i = 1:15
            unk = tophatfilter(unks[i], ff, 1.0 / dose(unks[i]))
            push!(res, filterfit(unk, fds))
        end

        # Compare against DTSA-II values
        @test isapprox(mean(values(res, oroi[1])), 0.4923, rtol = 0.003)
        @test isapprox(mean(values(res, siroi[1])), 0.0214, atol = 0.013)
        @test isapprox(mean(values(res, alroi[1])), 0.0281, atol = 0.001)
        @test isapprox(mean(values(res, caroi[1])), 0.1211, rtol = 0.0025)
        @test isapprox(mean(values(res, znroi[1])), 0.0700, rtol = 0.05)
        @test isapprox(mean(values(res, znroi[2])), 0.1115, atol = 0.0005)
        @test isapprox(mean(values(res, znroi[3])), 0.1231, rtol = 0.01)
        @test isapprox(mean(values(res, tiroi[1])), 0.0404, atol = 0.001)
        @test isapprox(mean(values(res, tiroi[2])), 0.06398, rtol = 0.0002)
        @test isapprox(mean(values(res, tiroi[3])), 0.064, rtol = 0.06)
        @test isapprox(mean(values(res, feroi[1])), 0.0, atol = 0.001)
        @test isapprox(mean(values(res, feroi[2])), 0.0, atol = 0.0004)
        @test isapprox(mean(values(res, feroi[3])), 0.0, atol = 0.001)
        @test isapprox(mean(values(res, geroi[1])), 0.1789, rtol = 0.01)
        @test isapprox(mean(values(res, geroi[2])), 0.2628, atol = 0.001)
        @test isapprox(mean(values(res, geroi[3])), 0.279, atol = 0.011)
    end
       @testset "ADM6005a - Refs" begin
        path = joinpath(@__DIR__, "ADM6005a spectra")
        unks = loadspectrum.(joinpath(path, "ADM-6005a_$(i).msa") for i = 1:15)
        det = matching(
            unks[1],
            128.0,
            110,
            Dict(
                Shell(1) => n"Be",
                Shell(2) => n"Sc",
                Shell(3) => n"Cs",
                Shell(4) => n"Pu",
            ),
        )
        ffp = references(
            [
                reference(n"Al", joinpath(path, "Al std.msa"), mat"Al"),
                reference(n"Ca", joinpath(path, "CaF2 std.msa"), mat"CaF2"),
                reference(n"Fe", joinpath(path, "Fe std.msa"), mat"Fe"),
                reference(n"Ge", joinpath(path, "Ge std.msa"), mat"Ge"),
                reference(n"Si", joinpath(path, "Si std.msa"), mat"Si"),
                reference(n"O", joinpath(path, "SiO2 std.msa"), mat"SiO2"),
                reference(n"Ti", joinpath(path, "Ti trimmed.msa"), mat"Ti"),
                reference(n"Zn", joinpath(path, "Zn std.msa"), mat"Zn"),
            ],
            det,
        )
        res = fit_spectra(unks, ffp)
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Al K-L3"))),
            0.0279,
            atol = 0.0001,
        )
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Ti K-L3"))),
            0.06399,
            atol = 0.0001,
        )
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Ge K-M3"))),
            0.274786,
            atol = 0.0001,
        )
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Zn K-M3"))),
            0.12103,
            atol = 0.0001,
        )
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Fe L3-M5"))),
            0.000290,
            atol = 0.00001,
        )
        @test isapprox(
            mean(values(res, findlabel(res[1], n"Fe K-L3"))),
            0.0003026,
            atol = 0.00001,
        )
    end

    # Check that the covariance of the filtered spectrum is calculated correctly as F*diagm(S)*transpose(F)
    @testset "Filtered covariance" begin
        spec = loadspectrum(joinpath(@__DIR__, "ADM6005a spectra", "ADM-6005a_1.msa"))
        det = matching(
            spec,
            128.0,
            110,
            Dict(
                Shell(1) => n"Be",
                Shell(2) => n"Sc",
                Shell(3) => n"Cs",
                Shell(4) => n"Pu",
            ),
        )
        filt = buildfilter(VariableWidthFilter, det)
        specdata = counts(spec)
        cov1 = [
            NeXLSpectrum.filteredcovar(filt, specdata, r, c) for
            r in eachindex(specdata), c in eachindex(specdata)
        ]
        filtd = NeXLSpectrum.filterdata(filt)
        cov2 = filtd * diagm(specdata) * transpose(filtd)
        # @show findmax(ii->abs(cov1[ii]-cov2[ii]), eachindex(cov1))
        @test all(
            isapprox(cov1[ii], cov2[ii], rtol = 1.0e-6, atol = 1.0e-12) for
            ii in eachindex(cov1)
        )
    end

    @testset "Repeated refs" begin
        path = joinpath(@__DIR__, "K412 spectra")
        fe = mat"Fe"
        efs = references(
                   [
                       reference(n"Ca", joinpath(path,"III-E K412[0][4].msa"), srm470_k412),
                       reference(n"Fe", joinpath(path,"III-E K412[0][4].msa"), srm470_k412),
                       reference(n"O", joinpath(path, "SiO2 std.msa"), mat"SiO2"),
                       reference(n"Al", joinpath(path, "Al2O3 std.msa"), mat"Al2O3"),
                       reference(n"Ca", joinpath(path, "CaF2 std.msa"), mat"CaF2"),
                       reference(n"Fe", joinpath(path, "Fe std.msa"), fe),
                       reference(n"Mg", joinpath(path, "MgO std.msa"), mat"MgO"),
                       reference(n"Si", joinpath(path, "SiO2 std.msa"), mat"SiO2"),
                   ],
                   132.0,
               )
        @test properties(efs.references[findfirst(r->n"Fe K-L3" in r.label.xrays, efs.references)].label)[:Composition] === srm470_k412
        @test properties(efs.references[findfirst(r->n"Fe K-M3" in r.label.xrays, efs.references)].label)[:Composition] === srm470_k412
        @test properties(efs.references[findfirst(r->n"Fe L3-M5" in r.label.xrays, efs.references)].label)[:Composition] === fe
        @test properties(efs.references[findfirst(r->n"Ca K-L3" in r.label.xrays, efs.references)].label)[:Composition] === srm470_k412
    end

    @testset "Warnings" begin
        s = loadspectrum(joinpath(@__DIR__, "Other", "K411 simulated.msa"))
        @test_logs ( :warn, "The spectrum \"Noisy[MC simulation of bulk K411] #1\" cannot be used as a reference for the ROI \"O K-L3 + 1 other\" due to 2 peak interferences.") 
            charXRayLabels(s, n"O", Set( ( n"C", n"O",n"Mg",n"Si",n"Ca",n"Fe")), simpleEDS(2048,10.0,0.0,132.0), 1.0e6, ampl=1.0e-5)
        @test_logs ( :warn, "The spectrum \"Noisy[MC simulation of bulk K411] #1\" cannot be used as a reference for the ROI \"Fe L3-M5 + 11 others\" due to 1 peak interference.") 
            charXRayLabels(s, n"Fe", Set( ( n"C", n"O",n"Mg",n"Si",n"Ca",n"Fe")), simpleEDS(2048,10.0,0.0,132.0), 1.0e6, ampl=1.0e-5)
    end

    @testset "Example 2" begin
        path = joinpath(@__DIR__, "Example 2")
        refs = references( [
            reference( [ n"Mg", n"Si", n"Ca", n"Fe" ], joinpath(path, "K411 std.msa"), srm470_k411)...,
            reference( n"O", joinpath(path,"MgO std.msa"), mat"MgO" ),
            reference( n"Fe", joinpath(path,"Fe std.msa"), mat"Fe" ),
            reference( n"Al", joinpath(path,"Al std.msa"), mat"Al" )
        ], 135.0)
        unk = loadspectrum(joinpath(path, "K412 unk.msa"))
        fr = fit_spectrum(unk, refs)
        qr = quantify(fr)
        @test isapprox(value(qr.comp[n"Al"]), 0.050776, atol=0.0001) 
        @test isapprox(value(qr.comp[n"Fe"]), 0.0783, atol=0.0001) 
        @test isapprox(value(qr.comp[n"Mg"]), 0.117958, atol=0.0001) 
        @test isapprox(value(qr.comp[n"O"]), 0.446058, atol=0.0001)

        df = asa(DataFrame, [ fr, fr ],  charOnly = false, withUnc = true, format = :normal)
        @test startswith(repr(df[1,:Spectra]),"\"K412-0[Mon Oct 17 16:11:17 2011]")
        @test ncol(df)==15
        @test nrow(df)==2
        @test isapprox(df[2,2],0.712819,atol=0.0001)
        @test isapprox(df[2,3],0.001413,atol=0.0001)

        df = asa(DataFrame, [ fr, ],  charOnly = false, withUnc = true, format = :pivot) 
        @test ncol(df)==3
        @test nrow(df)==7
        @test repr(df[1,:ROI])=="k[O K-L3 + 1 other, MgO]"
        @test isapprox(df[2,2], 0.0511784, atol=0.00001)
        @test isapprox(df[3,3], 0.00322127, atol=0.00001)

        df = asa(DataFrame, [ fr, fr ],  charOnly = false, withUnc = true, format = :long)
        @test ncol(df)==4 
        @test nrow(df)==14
        @test df[2,:ROI]=="k[Fe L3-M5 + 13 others, Fe]"
        @test isapprox(df[3,3], 1.3800169, atol=0.00001)
        @test isapprox(df[3,4],0.00322128, atol=0.00001)

        df = asa(DataFrame, fr, charOnly = false, material = srm470_k412, columns = ( :roi, :peakback, :counts, :dose))
        @test all(r->startswith(repr(r[:Spectrum]),"K412-0[Mon Oct 17 16:11:17 2011]"), eachrow(df))
        @test all(r->r[:LiveTime]==60.0,eachrow(df))
        @test all(r->r[:ProbeCurrent]==1.1978,eachrow(df))
        @test all(r->isapprox(r[:RealTime],69.97325,atol=0.0001),eachrow(df))
        @test df[1,:Start]==128
        @test df[1,:Stop]==173
        @test isapprox(df[1,:K], 0.0329847, atol=0.00001)
        @test isapprox(df[1,:dK], 0.0001588, atol=0.00001)
        @test isapprox(df[1,:Counts], 105385.0, atol=10.0)
        @test isapprox(df[1,:Back], 74516.9, atol=10.0)
        @test isapprox(df[1,:PtoB], 108.774, atol=0.001)
        @test isapprox(df[1,:KCalc], 0.032160115, atol=0.00001)
        @test isapprox(df[1,:KoKcalc], 1.0256414, atol=0.00002)
        @test isapprox(df[1,:RefCountsPernAs], 44455.97, atol=0.1)
        @test isapprox(df[1,:CountsPernAs], 1466.36, atol=0.1)
    end
    @testset "Example 2 - 32-bit" begin
        path = joinpath(@__DIR__, "Example 2")
        refs = references( [
            reference( [ n"Mg", n"Si", n"Ca", n"Fe" ], joinpath(path, "K411 std.msa"), srm470_k411)...,
            reference( n"O", joinpath(path,"MgO std.msa"), mat"MgO" ),
            reference( n"Fe", joinpath(path,"Fe std.msa"), mat"Fe" ),
            reference( n"Al", joinpath(path,"Al std.msa"), mat"Al" )
        ], 135.0, ftype=Float32) # This line is the only difference from "Example 2"
        unk = loadspectrum(joinpath(path, "K412 unk.msa"))
        fr = fit_spectrum(unk, refs)
        qr = quantify(fr)
        @test isapprox(value(qr.comp[n"Al"]), 0.0507768, atol=0.0001) 
        @test isapprox(value(qr.comp[n"Fe"]), 0.0783, atol=0.0001) 
        @test isapprox(value(qr.comp[n"Mg"]), 0.1179591, atol=0.0001) 
        @test isapprox(value(qr.comp[n"O"]), 0.4460586, atol=0.0001)

        df = asa(DataFrame, [ fr, fr ],  charOnly = false, withUnc = true, format = :normal)
        @test startswith(repr(df[1,:Spectra]),"\"K412-0[Mon Oct 17 16:11:17 2011]")
        @test ncol(df)==15
        @test nrow(df)==2
        @test isapprox(df[2,2],0.7128191,atol=0.0001)
        @test isapprox(df[2,3],0.001413,atol=0.0001)

        df = asa(DataFrame, [ fr, ],  charOnly = false, withUnc = true, format = :pivot) 
        @test ncol(df)==3 
        @test nrow(df)==7
        @test repr(df[1,:ROI])=="k[O K-L3 + 1 other, MgO]"
        @test isapprox(df[2,2], 0.051178406, atol=0.00001)
        @test isapprox(df[3,3], 0.003221281, atol=0.00001)

        df = asa(DataFrame, [ fr, fr ],  charOnly = false, withUnc = true, format = :long)
        @test ncol(df)==4 
        @test nrow(df)==14
        @test df[2,:ROI]=="k[Fe L3-M5 + 13 others, Fe]"
        @test isapprox(df[3,3], 1.3800191, atol=0.00001)
        @test isapprox(df[3,4], 0.0032212, atol=0.00001)

        df = asa(DataFrame, fr, charOnly = false, material = srm470_k412, columns = ( :roi, :peakback, :counts, :dose))
        @test all(r->startswith(repr(r[:Spectrum]),"K412-0[Mon Oct 17 16:11:17 2011]"), eachrow(df))
        @test all(r->r[:LiveTime]==60.0,eachrow(df))
        @test all(r->r[:ProbeCurrent]==1.1978,eachrow(df))
        @test all(r->isapprox(r[:RealTime], 69.97325,atol=0.0001),eachrow(df))
        @test df[1,:Start]==128
        @test df[1,:Stop]==173
        @test isapprox(df[1,:K], 0.03298476, atol=0.00001)
        @test isapprox(df[1,:dK], 0.0001588, atol=0.00001)
        @test isapprox(df[1,:Counts], 105385.09, atol=10.0)
        @test isapprox(df[1,:Back], 74516.9, atol=10.0)
        @test isapprox(df[1,:PtoB], 108.774, atol=0.001)
        @test isapprox(df[1,:KCalc], 0.032160, atol=0.00001)
        @test isapprox(df[1,:KoKcalc], 1.025642, atol=0.00002)
        @test isapprox(df[1,:RefCountsPernAs], 44455.9, atol=0.1)
        @test isapprox(df[1,:CountsPernAs], 1466.37, atol=0.1)
    end
end
