# X-ray Window

abstract type AbstractWindow end

abstract type WindowType end

struct MoxtekAP33 <: WindowType end
struct MoxtekAP5 <: WindowType end
struct AmetekC1 <: WindowType end
struct AmetekC2 <: WindowType end
struct BerylliumWindow <: WindowType end


Base.show(io::IO, wnd::AbstractWindow) = print(io, wnd.name)

struct NoWindow <: AbstractWindow
    name::String
    NoWindow() = new("No window")
end

NeXLCore.transmission(wnd::NoWindow, energy::Float64, angle::Float64 = π / 2) = 1.0


"""
    ModeledWindow(::Type{AmetekC1})
    ModeledWindow(::Type{AmetekC2})

Create modeled windows for the Ametek C1 or C2 Si₃N₄ windows according to the specifications here:
https://www.amptek.com/products/accessories-for-xrf-eds/c-series-low-energy-x-ray-windows#Specifications

    ModeledWindow(::Type{MoxtekAP33})
    ModeledWindow(::Type{MoxtekAP5})

Create modeled windows for Moxtek polymer windows.

    ModeledWindow(::Type{BerylliumWindow}; thickness = 5.0e-4)

Create a modeled window for a Be window (nominally 5 μm thick)
"""
struct ModeledWindow <: AbstractWindow
    name::String
    layers::Vector{Film}
    support::Film
    openfraction::Float64
end

NeXLCore.transmission(wnd::ModeledWindow, energy::Float64, angle::Float64 = π / 2) =
    (
        length(wnd.layers) > 0 ?
        mapreduce(lyr -> NeXLCore.transmission(lyr, energy, angle), *, wnd.layers) : 1.0
    ) * (
        wnd.openfraction +
        (1.0 - wnd.openfraction) * NeXLCore.transmission(wnd.support, energy, angle)
    )

function ModeledWindow(::Type{MoxtekAP33})
    support, openarea = Film(pure(n"Si"), 0.038), 0.77
    paralene = Film(parse(Material, "C10H8O4N", density = 1.39), 3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    ModeledWindow("Moxtek AP3.3 model", [aluminum, paralene], support, openarea)
end

function ModeledWindow(::Type{MoxtekAP5})
    support, openarea = Film(pure(n"C"), 0.0265), 0.78
    paralene = Film(parse(Material, "C10H8O4N", density = 1.39), 3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    ModeledWindow("Moxtek AP5 model", [aluminum, paralene], support, openarea)
end

function ModeledWindow(::Type{BerylliumWindow}; thickness = 5.0e-4)
    support, openarea = Film(pure(n"Be"), thickness), 0.0
    ModeledWindow("$(thickness*1.0e4) μm Be window", [], support, openarea)
end

function ModeledWindow(::Type{AmetekC1})
    si3n4 = Film(parse(Material, "Si3N4", density = 3.44), 150.0e-7) # 150 nm of Si₃N₄
    al = Film(pure(n"Al"), 250.0e-7) # 250 nm of Al
    support, openarea = Film(pure(n"Si"), 1.50e-4), 0.80 # 80% open area, 15 μm Si thickness support
    ModeledWindow("AMETEK C1", [si3n4, al], support, openarea)
end

function ModeledWindow(::Type{AmetekC2})
    si3n4 = Film(parse(Material, "Si3N4", density = 3.44), 40.0e-7) # 40 nm of Si₃N₄
    al = Film(pure(n"Al"), 30.0e-7) # 30 nm of Al
    support, openarea = Film(pure(n"Si"), 1.50e-4), 0.80 # 15 μm thickness Si support
    ModeledWindow("AMETEK C2", [si3n4, al], support, openarea)
end

"""
    TabulatedWindow(::Type{MoxtekAP33})
    TabulatedWindow(::Type{MoxtekAP5})
    TabulatedWindow(::Type{AmetekC1})
    TabulatedWindow(::Type{AmetekC2})

Construct tabulated window models for the Moxtek AP3.3 and AP5 windows and Ametek C1 and C2 windows.
"""
struct TabulatedWindow <: AbstractWindow
    name::String
    interpolation::AbstractInterpolation
    extrapolation::ModeledWindow
    match::Float64
end

function NeXLCore.transmission(
    wnd::TabulatedWindow,
    energy::Float64,
    angle::Float64 = π / 2,
)
    bnds = bounds(parent(wnd.interpolation))[1]
    return (energy >= bnds[1]) && (energy < bnds[2]) ? wnd.interpolation(energy) : #
           wnd.match * transmission(wnd.extrapolation, energy, angle)
end

function TabulatedWindow(::Type{MoxtekAP33})
    data = CSV.read(
        joinpath(@__DIR__, "data", "AP3_3_mod.csv"),
        DataFrame,
        comment = "//"
    )
    inter = LinearInterpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(MoxtekAP33)
    TabulatedWindow("Moxtek AP3.3", inter, extra, 1.01582935271)
end

function TabulatedWindow(::Type{MoxtekAP5})
    data = CSV.read(joinpath(@__DIR__, "data", "AP5.csv"), DataFrame)
    inter = LinearInterpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(MoxtekAP5)
    TabulatedWindow("Moxtek AP5", inter, extra, 1.0)
end

function TabulatedWindow(::Type{AmetekC1})
    data = CSV.read(joinpath(@__DIR__, "data", "AMETEK Si3N4 C1.csv"), DataFrame, header=2)
    inter = LinearInterpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(AmetekC1)
    TabulatedWindow("AMETEK C1 Si₃N₄", inter, extra, 1.0)
end

function TabulatedWindow(::Type{AmetekC2})
    data = CSV.read(joinpath(@__DIR__, "data", "AMETEK Si3N4 C2.csv"), DataFrame, header=2)
    inter = LinearInterpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(AmetekC2)
    TabulatedWindow("AMETEK C2 Si₃N₄", inter, extra, 1.0)
end

