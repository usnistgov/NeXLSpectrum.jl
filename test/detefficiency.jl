using Test

@testset "Efficiency" begin
    sdd = SDDEfficiency(ModeledWindow(MoxtekAP33()))
    sili = SiLiEfficiency(ModeledWindow(BerylliumWindow(10.0e-4)))

    @test isapprox(efficiency(sdd, 300.0), 0.152, atol = 0.001)
    @test isapprox(efficiency(sdd, 3000.0), 0.753, atol = 0.001)
    @test isapprox(efficiency(sdd, 30000.0), 0.091, atol = 0.001)

    @test isapprox(efficiency(sili, 300.0), 0.000, atol = 0.001)
    @test isapprox(efficiency(sili, 3000.0), 0.957, atol = 0.001)
    @test isapprox(efficiency(sili, 30000.0), 0.485, atol = 0.001)
end

@testset "Simulate" begin
    basepath = joinpath(@__DIR__, "ADM6005a spectra")
    sp = loadspectrum(joinpath(basepath,"ADM-6005a_1.msa"))
    det = matching(sp, 130.0)
    refs = references( [
        reference(n"O", joinpath(basepath,"SiO2 std.msa"), mat"SiO2"),
        reference(n"Si", joinpath(basepath,"SiO2 std.msa"), mat"SiO2"),
        reference(n"Al", joinpath(basepath,"Al std.msa"), mat"Al"),
        reference(n"Ca", joinpath(basepath,"CaF2 std.msa"), mat"CaF2"),
        reference(n"Ge", joinpath(basepath,"Ge std.msa"), mat"Ge"),
        reference(n"Ti", joinpath(basepath,"Ti trimmed.msa"), mat"Ti"),
        reference(n"Zn", joinpath(basepath,"Zn std.msa"), mat"Zn"),
    ], det)
    k = fit_spectrum(sp, refs)
    q1, q2 = quantify(k), quantify(k, kro = SimpleKRatioOptimizer(1.5, [ n"Ge L3-M5" ]))


    eff = SDDEfficiency(ModeledWindow(MoxtekAP33()))
    resp = detectorresponse(det, eff)
    sim1 = NeXLSpectrum.simulate(nonneg(q1.comp), dose(sp), sp[:BeamEnergy], sp[:TakeOffAngle], 40.0/(72.0^2), det, resp)
    sim2 = NeXLSpectrum.simulate(nonneg(q1.comp), dose(sp), sp[:BeamEnergy], sp[:TakeOffAngle], 40.0/(72.0^2), det, resp)
end