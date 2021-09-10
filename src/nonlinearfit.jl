# Implements a non-linear fit

using Polynomials
using NeXLUncertainties: LabeledValues, Label
using NeXLSpectrum
using Gadfly

abstract type FitFunction end

struct CoeffN <: Label
    n::Int
end 
struct FWHMLabel <: Label end
struct PeakI <: Label
    xray::CharXRay
end

struct PeakShapeFF <: FitFunction
    material::Material
    order::Int

    PeakShapeFF(mat::Material, cxrs::AbstractVector{CharXRay}, order::Int) = new(mat, order)
end

function fitfunc(ps::PeakShapeFF, args::LabeledValues)::Function
    det = BasicEDS(
        length(chs), 
        channelcount::Int,
        LinearEnergyScale(args[CoeffN(0)], args[CoeffN(1)]),
        MnKaResolution(args[FWHMLabel()]),
        10
    )
    
    resp = detectorresponse(det, SDDEfficiency(AP33Model()))

    @show mnka
    imaxs = map(ps.xrays) do xray
        args[PeakI(xray)]*weight(xray)
    end
    @show imaxs
    exs = map(xray->energy(xray), ps.xrays)
    @show exs
    return ch -> begin
        e = efunc(ch)
        sum(i->imaxs[i]*profile(e, exs[i], mnka), eachindex(imaxs))
    end
end


function plotit()
    chs = 1:2048


    psff = PeakShapeFF(characteristic(n"Fe",alltransitions), 1)

    cxrs = characteristic(n"Fe", alltransitions)
    args = LabeledValues(
        map(cxr->PeakI(cxr)=> 1.0, cxrs)...,
        CoeffN(0) => 0.0,
        CoeffN(1) => 10.0,
        FWHMLabel() =>132.0
        )

    ff = fitfunc(psff, args)
    plot(x=chs, y = ff.(chs), Geom.line)    
end