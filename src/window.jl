"""
`AbstractWindow` is the base type for `TabulatedWindow`, `ModeledWindow` and `NoWindow`.
You can view the window transmission function using `Gadfly` and the method:

    plot(wind::Union{AbstractWindow, AbstractArray{<:AbstractWindow}}; xmax=20.0e3, angle=π/2, style=NeXLSpectrumStyle)

Window types are identified by types based on `WindowType` like `MoxtekAP33`, `MoxtekAP5`, `AmptekC1`, `AmptekC2`, `BerylliumWindow`.

An implementation of an `AbstractWindow` is instantiate using code like `TabulatedWindow(MoxtekAP33())` or `ModeledWindow(AmptekC1())`.
In addition, there is a `NoWindow()` type that implements a 100% transparent window.

`AbstractWindow` types primarily implement `NeXLCore.transmission(wnd::AbstractWindow, energy::Float64, angle::Float64 = π / 2)`
"""
abstract type AbstractWindow end
Base.show(io::IO, wnd::AbstractWindow) = print(io, name(wnd))

"""
    NoWindow

A 100% transparent window.
"""
struct NoWindow <: AbstractWindow end
NeXLCore.name(::NoWindow) = "No window"
NeXLCore.transmission(wnd::NoWindow, energy::Float64, angle::Float64=π / 2) = 1.0


"""
`WindowType` is the abstract type for models or types of X-ray windows.  These types distinguish between
the window type and the implementation.  Often, there are both tabulated transmission functions from
the vendor and calculated transmission functions based on the construction of the window.
Use the `TabulatedWindow` or `ModeledWindow` to instantiate `AbstractWindow` which implements the
`transmission(wnd::AbstractWindow, energy::Float64, angle::Float64 = π / 2)` method.

Predefined `WindowType`s are `MoxtekAP33`, `MoxtekAP5`, `AmptekC1`, `AmptekC2`, `BerylliumWindow`.
"""
abstract type WindowType end
struct MoxtekAP33 <: WindowType end
NeXLCore.name(::MoxtekAP33) = "Moxtek AP3.3"
struct MoxtekAP5 <: WindowType end
NeXLCore.name(::MoxtekAP5) = "Moxtek AP5"

"""
Create modeled windows for the Ametek C1 Si₃N₄ windows according to the specifications here:
    https://www.amptek.com/products/accessories-for-xrf-eds/c-series-low-energy-x-ray-windows#Specifications   

Light-tight (solar-blind) window. Often used for environmental XRF units.
"""
struct AmptekC1 <: WindowType end
NeXLCore.name(::AmptekC1) = "AMPTEK C1 Si₃N₄"

"""
Create modeled windows for the Ametek C2 Si₃N₄ windows according to the specifications here:
    https://www.amptek.com/products/accessories-for-xrf-eds/c-series-low-energy-x-ray-windows#Specifications   

Light transparent window often used for SEM detectors.
"""
struct AmptekC2 <: WindowType end
NeXLCore.name(::AmptekC2) = "AMPTEK C2 Si₃N₄"

# Aliases: Amptek was bought by Ametek (confusing???)
const AmetekC1 = AmptekC1;
const AmetekC2 = AmptekC2;


"""
    BerylliumWindow(thickness)

Create a window of pure Be of the specified thickness (in cm)
"""
struct BerylliumWindow <: WindowType
    thickness::Float64

    BerylliumWindow(thickness=5.0e-4) = new(thickness)
end
NeXLCore.name(bw::BerylliumWindow) = "$(bw.thickness*1.0e4) μm Beryllium"


"""
    ModeledWindow(wt::WindowType)

This type models a window using the materials and thicknesses provided by the vendor.
If accomodates grids by a simplistic mechanism of assuming an open area.
"""
struct ModeledWindow <: AbstractWindow
    type::WindowType
    layers::Vector{Film}
    support::Film
    openfraction::Float64
end
NeXLCore.name(mw::ModeledWindow) = "$(name(mw.type)) - Modeled"

function NeXLCore.transmission(
    wnd::ModeledWindow,
    energy::Float64,
    angle::Float64=π / 2
)
    lt = mapreduce(lyr -> NeXLCore.transmission(lyr, energy, angle), *, wnd.layers, init=1.0)
    lt * (wnd.openfraction + (1.0 - wnd.openfraction) * NeXLCore.transmission(wnd.support, energy, angle))
end
NeXLCore.transmission(wnd::AbstractWindow, cxr::CharXRay, angle::Float64=π / 2) = transmission(wnd, energy(cxr), angle)


function ModeledWindow(wt::MoxtekAP33)
    support, openarea = Film(pure(n"Si"), 0.038), 0.77
    paralene = Film(parse(Material, "C10H8O4N", density=1.39), 3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    ModeledWindow(wt, [aluminum, paralene], support, openarea)
end

function ModeledWindow(wt::MoxtekAP5)
    support, openarea = Film(pure(n"C"), 0.0265), 0.78
    paralene = Film(parse(Material, "C10H8O4N", density=1.39), 3.0e-5)
    aluminum = Film(pure(n"Al"), 4.0e-6)
    ModeledWindow(wt, [aluminum, paralene], support, openarea)
end

function ModeledWindow(wt::BerylliumWindow)
    support, openarea = Film(pure(n"Be"), wt.thickness), 0.0
    ModeledWindow(wt, [], support, openarea)
end

function ModeledWindow(wt::AmptekC1, openarea=0.80)
    si3n4 = Film(parse(Material, "Si3N4", density=3.44), 150.0e-7) # 150 nm of Si₃N₄
    al = Film(pure(n"Al"), 250.0e-7) # 250 nm of Al
    support = Film(pure(n"Si"), 1.50e-3) # 80% open area, 15 μm Si thickness support
    ModeledWindow(wt, [si3n4, al], support, openarea)
end

function ModeledWindow(wt::AmptekC2, openarea=0.80)
    si3n4 = Film(parse(Material, "Si3N4", density=3.44), 40.0e-7) # 40 nm of Si₃N₄
    al = Film(pure(n"Al"), 30.0e-7) # 30 nm of Al
    support = Film(pure(n"Si"), 1.50e-3) # 15 μm thickness Si support
    ModeledWindow(wt, [si3n4, al], support, openarea)
end

"""
    TabulatedWindow(wt::WindowType)

Construct a model of the window transmission based on vendor-supplied tabulations of transparency. At the
end of the user supplied data, the transmission function is extended by matching the ModeledWindow
with the tabulation.
"""
struct TabulatedWindow <: AbstractWindow
    type::WindowType
    interpolation::AbstractInterpolation
    extrapolation::ModeledWindow
    match::Float64 # Used to match tabulated and modeled transmission functions.
end
NeXLCore.name(mw::TabulatedWindow) = "$(name(mw.type)) - Tabulated"

function NeXLCore.transmission(
    wnd::TabulatedWindow,
    energy::Float64,
    angle::Float64=π / 2,
)
    bnds = bounds(parent(wnd.interpolation))[1]
    return (energy >= bnds[1]) && (energy < bnds[2]) ? wnd.interpolation(energy) : #
           wnd.match * transmission(wnd.extrapolation, energy, angle)
end

function TabulatedWindow(wt::MoxtekAP33)
    data = CSV.read(joinpath(@__DIR__, "data", "AP3_3_mod.csv"), DataFrame, header=3, comment="//")
    inter = linear_interpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(wt)
    match = inter(data[end, 1]) / transmission(extra, data[end, 1], π / 2)
    TabulatedWindow(wt, inter, extra, match)
end

function TabulatedWindow(wt::MoxtekAP5)
    data = CSV.read(joinpath(@__DIR__, "data", "AP5.csv"), DataFrame)
    inter = linear_interpolation(data[:, 1], data[:, 2])
    extra = ModeledWindow(wt)
    match = inter(data[end, 1]) / transmission(extra, data[end, 1], π / 2)
    TabulatedWindow(wt, inter, extra, match)
end

function TabulatedWindow(wt::AmptekC1)
    data = CSV.read(joinpath(@__DIR__, "data", "AMETEK Si3N4 C1.csv"), DataFrame, header=2)
    inter = linear_interpolation(1000.0 * data[:, 1], data[:, 2])
    extra = ModeledWindow(wt)
    match = inter(1000.0 * data[end, 1]) / transmission(extra, 1000.0 * data[end, 1], π / 2)
    TabulatedWindow(wt, inter, extra, match)
end

function TabulatedWindow(wt::AmptekC2)
    data = CSV.read(joinpath(@__DIR__, "data", "AMETEK Si3N4 C2.csv"), DataFrame, header=2)
    inter = linear_interpolation(1000.0 * data[:, 1], data[:, 2])
    extra = ModeledWindow(wt)
    match = inter(1000.0 * data[end, 1]) / transmission(extra, 1000.0 * data[end, 1], π / 2)
    TabulatedWindow(wt, inter, extra, match)
end
