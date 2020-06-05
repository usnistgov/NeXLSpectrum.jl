# X-ray Window
using NeXLCore
using CSV
using Interpolations: LinearInterpolation, AbstractInterpolation, bounds

abstract type AbstractWindow end

Base.show(io::IO, wnd::AbstractWindow) = print(io, wnd.name)

struct NoWindow <: AbstractWindow
    name::String
    NoWindow() = new("No window")
end

NeXLCore.transmission(wnd::NoWindow, energy::Float64, angle::Float64=π/2) = 1.0

struct LayerWindow <: AbstractWindow
    name::String
    layers::Vector{Film}
    support::Film
    openfraction::Float64
end

NeXLCore.transmission(wnd::LayerWindow, energy::Float64, angle::Float64=π/2) =
    (length(wnd.layers) > 0 ? mapreduce(lyr->NeXLCore.transmission(lyr, energy, angle), *, wnd.layers) : 1.0) *
        (wnd.openfraction + (1.0-wnd.openfraction)*NeXLCore.transmission(wnd.support, energy, angle))

"""
    AP33Model()

Construct a modeled window for the Moxtek AP3.3 window.
"""

function AP33Model()
    support, openarea = Film(pure(n"Si"), 0.038), 0.77
    paralene = Film(parse(Material, "C10H8O4N", density=1.39),3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    return LayerWindow("Moxtek AP3.3 model", [aluminum, paralene], support, openarea)
end

function AP5Model()
    support, openarea = Film(pure(n"C"), 0.0265), 0.78
    paralene = Film(parse(Material, "C10H8O4N", density=1.39),3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    return LayerWindow("Moxtek AP5 model", [aluminum, paralene], support, openarea)
end

"""
    Beryllium(thickness=5.0e-4)

Construct a beryllium window.
"""
function Beryllium(thickness=5.0e-4)
    support, openarea = Film(pure(n"Be"), thickness), 0.0
    return LayerWindow("$(thickness*1.0e4) μm Be window", [ ], support, openarea)
end

"""
    function AmptekC1()
    function AmptekC2()

Create modeled windows for the Amptek C1 or C2 Si₃N₄ windows according to the specifications here:
https://www.amptek.com/products/accessories-for-xrf-eds/c-series-low-energy-x-ray-windows#Specifications
"""

function AmptekC1()
    si3n4 = Film(parse(Material, "Si3N4", density=3.44), 150.0e-7)
    al = Film(pure(n"Al"), 250.0e-7)
    support, openarea = Film(pure(n"Si"), 1.50e-4), 0.80
    return LayerWindow("Amptek C1", [ si3n4, al ], support, openarea)
end

function AmptekC2()
    si3n4 = Film(parse(Material, "Si3N4", density=3.44), 40.0e-7)
    al = Film(pure(n"Al"), 30.0e-7)
    support, openarea = Film(pure(n"Si"), 1.50e-4), 0.80
    return LayerWindow("Amptek C2", [ si3n4, al ], support, openarea)
end

struct TabulatedWindow <: AbstractWindow
    name::String
    interpolation::AbstractInterpolation
    extrapolation::LayerWindow
    match::Float64
end

function NeXLCore.transmission(wnd::TabulatedWindow, energy::Float64, angle::Float64=π/2)
    bnds = bounds(parent(wnd.interpolation))[1]
    return (energy >= bnds[1]) && (energy < bnds[2]) ? wnd.interpolation(energy) : #
        wnd.match*transmission(wnd.extrapolation,energy,angle)
end

"""
    AP33Tabulation()
    AP5Tabulation()

Construct tabulated window models for the Moxtek AP3.3 and AP5 windows.
"""
function AP33Tabulation()
    data = CSV.read(joinpath(dirname(pathof(@__MODULE__)), "AP3_3_mod.csv"), skipto=3)
    inter = LinearInterpolation(data[:,1], data[:,2])
    extra = AP33Model()
    return TabulatedWindow("Moxtek AP3.3", inter, extra, 1.01582935271)
end


function AP5Tabulation()
    data = CSV.read(joinpath(dirname(pathof(@__MODULE__)), "AP5.csv"))
    inter = LinearInterpolation(data[:,1], data[:,2])
    extra = AP5Model()
    return TabulatedWindow("Moxtek AP5",inter, extra, 1.0)
end
""
