using NeXLSpectrum
using Test

@testset "XRF Spectrum" begin
    fn1 = joinpath("Other","20170510_D1_Mo_50kv_100mA_acryl 0.spx")
    sp1a = loadspectrum(fn1)
    #sp1b = FileIO.load(fn1)
    @test sp1a[:RealTime] == 310.610
    @test sp1a[:LiveTime] == 300.294
    @test sp1a[:BrukerThroughput] == 130000
    @test sp1a[:DetectorModel] == "Bruker XFlash 430"
    @test sp1a[:DetectorSerialNumber] == "11881_0239"
    @test sp1a[:DetectorThickness] == 0.1*0.45 # cm
    @test sp1a[:DeadLayerThickness] == 0.1*0.029 # ???
    @test sp1a[:BeamEnergy] == 1000.0 * 50.0
    @test sp1a[:XRFTubeAnode] == n"Mo"
    @test sp1a[:ProbeCurrent] == 1000.0*99
    @test sp1a[:XRFTubeIncidentAngle] == deg2rad(84.0)
    @test sp1a[:XRFTubeTakeOffAngle] == deg2rad(6.0)
    @test sp1a[:XRFExcitationAngle] == deg2rad(50.0)
    @test sp1a[:Elevation] == deg2rad(50.0)
    @test sp1a[:TakeOffAngle] == deg2rad(50.0)
    @test sp1a[:XRFExcitationPathLength] == 0.1*10.0
    @test sp1a[:XRFDetectionPathLength] == 0.1*20.0
    @test sp1a[:DetectorSolidAngle] == 0.0065
    @test sp1a[:ChamberPressure] == 20.0
    @test sp1a[:ChamberAtmosphere] == "Air"
    @test sp1a[:XRFSampleTilt] == 0.0
    f1 = sp1a[:XRFTubeWindow]
    @test isequal(f1.material, pure(n"Be"))
    @test isequal(f1.thickness ,200.0e-4)
end

@testset "SEM/EDS" begin
    fn2 = joinpath("Other","Albite_ChMxd_20kV0p5nA130kHzTC_5p9kHzOCR_4ks.spx")
    sp2a = loadspectrum(fn2)
    #sp2b = FileIO.load(fn2)
    @test sp2a[:RealTime] == 4028.702
    @test sp2a[:LiveTime] == 4000.005
    @test sp2a[:BrukerThroughput] == 130000
    @test sp2a[:DetectorModel] == "Bruker XFlash 6|30"
    @test sp2a[:DetectorSerialNumber] == "13773"
    @test sp2a[:DetectorThickness] == 0.1*0.45
    @test sp2a[:DeadLayerThickness] == 0.1*0.029
    @test sp2a[:Elevation]==deg2rad(35.0)
    win = sp2a[:Window]
    @test length(win)==5
    @test isapprox(win[1], Film(pure(n"B"),1.3E1 * 1.0e-7))
    @test isapprox(win[2], Film(pure(n"C"),1.45E1 * 1.0e-7))
    @test isapprox(win[3], Film(pure(n"N"),4.5 * 1.0e-7))
    @test isapprox(win[4], Film(pure(n"O"),8.5 * 1.0e-7))
    @test isapprox(win[5], Film(pure(n"Al"),3.5 * 1.0e-7))
    @test sp2a[:BeamEnergy] == 20.0e3

    @test n"C" in sp2a[:Elements]
    @test n"O" in sp2a[:Elements]
    @test n"Na" in sp2a[:Elements]
    @test n"Al" in sp2a[:Elements]
    @test n"Si" in sp2a[:Elements]
    @test sp2a[:TakeOffAngle] == deg2rad(35.0)
    @test sp2a[:WorkingDistance] == 1.1
end
