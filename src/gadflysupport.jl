using .Gadfly
using Colors

const NeXLPalette = ( # This palette - https://flatuicolors.com/palette/nl
	RGB(234/255, 32/255, 39/255), RGB(27/255, 20/255, 100/255),
	RGB(0/255, 98/255, 102/255), RGB(87/255, 88/255, 187/255),
	RGB(11/2551, 30/255, 81/255), RGB(247/255, 159/255, 31/255),
	RGB(18/255, 137/255, 167/255), RGB(163/255, 203/255, 56/255),
	RGB(217/255, 128/255, 250/255), RGB(181/255, 52/255, 113/255),
	RGB(238/255, 90/255, 36/255), RGB(6/255, 82/255, 221/255),
	RGB(0/255, 148/255, 50/255), RGB(153/255, 128/255, 250/255),
	RGB(131/255, 52/255, 113/255), RGB(255/255, 195/255, 18/255),
	RGB(18/255, 203/255, 196/255), RGB(196/255, 229/255, 56/255),
	RGB(253/255, 167/255, 223/255), RGB(237/255, 76/255, 103/255) )

"""
    plot(spec::Spectrum)

Plot a Spectrum using Gadfly
"""
function Gadfly.plot(spec::Spectrum, klmLines=[]; xmin=0.0, xmax=nothing)::Plot
    maxCh = channel(get(spec, :BeamEnergy, energy(length(spec),spec)),spec)
    norm = dose(spec)
    data = (isequal(norm, missing) ? 1.0 : 1.0/norm) * spec.counts[1:maxCh]
    if isnothing(xmax)
        xmax=energy(maxCh,spec)
    end
    plot(x=map(ch->energy(ch, spec), 1:maxCh), y=data, Geom.step,
        Guide.XLabel("Energy (eV)"),Guide.YLabel(isequal(norm, missing) ? "Counts" : "Counts/(nA⋅s)"),
        Guide.Title(spec.properties[:Name]),
        Scale.x_continuous(format=:plain), Scale.y_continuous(format=:plain),
        Coord.Cartesian(ymin=0,ymax=1.05*maximum(data),xmin=xmin,xmax=xmax))
end

"""
    plot(specs::Vector{Spectrum})

Plot a multiple spectra on a single plot using Gadfly.
"""
function Gadfly.plot(specs::Vector{Spectrum}; xmin=0.0, xmax=nothing)::Plot
    maxI, maxE = 16, 1.0e3
    layers=[]
    for (i, spec) in enumerate(specs)
        chs = max(1,channel(xmin,spec)):channel(get(spec, :BeamEnergy, energy(length(spec),spec)),spec)
        maxI = maximum( [ maxI, maximum(spec.counts[chs]) ] )
        maxE = maximum( [ maxE, energy(chs.stop,spec) ] )
        push!(layers, layer(x=energyscale(spec)[chs], y=spec.counts[chs], Geom.step, Theme(default_color=NeXLPalette[(i-1) % length(NeXLPalette)+1] )))
    end
    plot(layers...,
        Guide.XLabel("Energy (eV)"), Guide.YLabel("Counts"),
        Scale.x_continuous(format=:plain), Scale.y_continuous(format=:plain),
        Coord.Cartesian(ymin=0, ymax=1.05*maxI, xmin=xmin, xmax= isnothing(xmax) ? maxE : xmax))
end

"""
    plot(fd::FilteredData)

Plot the fitting filter and the spectrum from which it was derived.  Vertical red
lines represent the principle ROC.
"""
function Gadfly.plot(fd::FilteredDatum)
    p1=layer(x=fd.ffroi, y=fd.filtered, Geom.step, Theme(default_color="darkcyan"),
        xintercept=[fd.roi.start, fd.roi.stop], Geom.vline(color=["firebrick","firebrick"]))
    p2=layer(x=fd.ffroi, y=fd.data, Geom.step, Theme(default_color="blueviolet"))
    plot(p1,p2,Coord.cartesian(xmin=fd.ffroi.start, xmax=fd.ffroi.stop),
            Guide.XLabel("Channels"), Guide.YLabel("Counts"), Guide.title(repr(fd)))
end
