# Model detector efficiency (Window +  Detector Crystal)

struct DetectorEfficiency
    name::String
    window::AbstractWindow
    surface::Vector{Film}
    active::Film
end

Base.show(io::IO, de::DetectorEfficiency) = print(io, "$(de.name)[$(de.window)]")

SDDEfficiency(window::AbstractWindow; thickness=0.0370, deadlayer=30.0e-7, entrance=Film(pure(n"Al"), 10.0e-7)) =
    DetectorEfficiency("SDD", window, [ entrance, Film(pure(n"Si"), deadlayer)], Film(pure(n"Si"), thickness))

SiLiEfficiency(window::AbstractWindow;thickness=0.250, deadlayer=30.0e-7, entrance=Film(pure(n"Al"), 10.0e-7)) =
    DetectorEfficiency("Si(Li)", window, [ entrance, Film(pure(n"Si"), deadlayer)], Film(pure(n"Si"), thickness))

efficiency(aa::DetectorEfficiency, energy::Float64, angle::Float64=π/2) =
    transmission(aa.window, energy, angle)*(1.0-transmission(aa.active, energy, angle))*mapreduce(lyr->transmission(lyr, energy, angle), *, aa.surface)
