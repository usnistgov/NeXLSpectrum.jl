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
    Gadfly.plot(spec::Spectrum; klms=[], xmin=0.0, xmax=nothing)::Plot

Plot a Spectrum using Gadfly.  klms is a Vector of CharXRays or Elements.
"""
Gadfly.plot(spec::Spectrum; klms=[], edges=[], autoklms=false, xmin=0.0, xmax=missing, norm=:None, lld=missing, yscale=1.05, ytransform=identity)::Plot =
	plot([spec], klms=klms, edges=edges, autoklms=autoklms, xmin=xmin, xmax=xmax, norm=norm, lld=lld, yscale=yscale, ytransform=identity)

"""
    plot(
	    specs::AbstractVector{Spectrum};
	    klms=[],
	    edges=[],
	    autoklms = false,
	    xmin=0.0,
	    xmax=missing,
	    norm=:None,
	    lld=missing,
	    yscale=1.05,
	    ytransform = identity
    )::Plot

Plot a multiple spectra on a single plot using Gadfly.

    norm = :None|:Sum|:Dose|:Peak|:DoseWidth
	klms = [ Element &| CharXRay ]
	edges = [ Element &| AtomicShell ]
	autoklms = false # Add KLMs based on elements in spectra
	lld = missing # missing defaults to 100.0  eV (low level discriminator for peak and intensity scaling)
	xmin = 0.0 # Min energy (eV)
	xmax = missing # Max energy (eV) (defaults to max(:BeamEnergy))
	yscale = 1.05 # Fraction of max intensity for ymax
	ytransform = identity | log10 | sqrt | ???
"""
function Gadfly.plot(
	specs::AbstractVector{Spectrum};
	klms=[],
	edges=[],
	autoklms = false,
	xmin=0.0,
	xmax=missing,
	norm=:None,
	lld=missing,
	yscale=1.05,
	ytransform = identity
)::Plot
	dlld(spec) = channel(ismissing(lld) ? (haskey(spec,:Detector) ? lld(spec[:Detector]) : 100.0) : lld, spec)
	# The normalize functions return a Vector{Vector{Float64}}
	normalizeDoseWidth(specs::AbstractVector{Spectrum}, def=1.0) =
		normalizedosewidth.(specs)
	normalizeDose(specs::AbstractVector{Spectrum}, def=1.0) =
		collect( (def/dose(sp))*counts(sp, Float64) for sp in specs)
	normalizeSum(specs::AbstractVector{Spectrum}, total=1.0e6) =
		collect( (total/integrate(sp,dlld(sp):length(sp))) * counts(sp, Float64) for sp in specs)
	normalizePeak(specs::AbstractVector{Spectrum}, height=100.0) =
		collect( (height/findmax(sp, dlld(sp):length(sp))[1]) * counts(sp, Float64) for sp in specs)
	normalizeNone(specs::AbstractVector{Spectrum}) =
		collect( counts(sp, Float64) for sp in specs)
	function klmLayer(specdata, cxrs::AbstractArray{CharXRay})
	    d=Dict{Any,Array{CharXRay}}()
	    for cxr in cxrs
	        d[(element(cxr),family(cxr))] = push!(get(d, (element(cxr), family(cxr)), []), cxr)
	    end
	    x, y, label = [], [], []
	    for cs in values(d)
	        br=brightest(cs)
			ich = maximum(specdata[i][channel(energy(br), spec)] for (i,spec) in enumerate(specs))
	        if ich > 0
	            for c in cs
	                push!(x, energy(c))
	                push!(y, ytransform(ich*weight(c)))
					push!(label, weight(c) > 0.1 ? "$(element(c).symbol)" : "")
	            end
	        end
	    end
	    return layer(x=x, y=y, label=label, Geom.hair, Geom.point, Geom.label(position=:above), Theme(default_color="gray" ))
	end
	function edgeLayer(maxI, ashs::AbstractArray{AtomicShell})
		maxCapacity(ashs) = largest
	    d=Dict{Any,Array{AtomicShell}}()
	    for ash in ashs
	        d[(element(ash),family(ash))] = push!(get(d, (element(ash), family(ash)), []), ash)
	    end
	    x, y, label = [], [], []
	    for ass in values(d)
	        br=ass[findmax(capacity.(ass))[2]]
            for ash in ass
                push!(x, energy(ash))
                push!(y, ytransform(maxI*capacity(ash)/capacity(br)))
                push!(label, "$(ash)")
            end
	    end
	    return layer(x=x, y=y, label=label, Geom.hair, Geom.label(position=:right), Theme(default_color="lightgray" ))
	end
	if norm==:Dose
		specdata=normalizeDose(specs)
		ylbl = "Counts/(nA⋅s)"
	elseif norm==:Sum
		specdata=normalizeSum(specs)
		ylbl = "Normalized (Σ=10⁶)"
	elseif norm==:Peak
		specdata=normalizePeak(specs)
		ylbl = "Peak (%)"
	elseif norm==:DoseWidth
		specdata=normalizeDoseWidth(specs)
		ylbl = "Counts/(nA⋅s⋅eV)"
	else
		specdata=normalizeNone(specs)
		ylbl = "Counts"
	end
    maxI, maxE = 16, 1.0e3
    names, layers, colors=[], [],[]
    for (i, spec) in enumerate(specs)
		mE = ismissing(xmax) ? get(spec, :BeamEnergy, energy(length(spec), spec)) : xmax
        chs = max(1,channel(xmin,spec)):channel(mE,spec)
		mchs = max(chs.start, dlld(spec)):chs.stop  # Ignore zero strobe...
        maxI = max(maxI, maximum(specdata[i][mchs]))
        maxE = max(maxE, mE)
		clr = NeXLSpectrum.NeXLPalette[(i-1) % length(NeXLSpectrum.NeXLPalette)+1]
		push!(names, spec[:Name])
		push!(colors, clr)
        push!(layers, layer(x=energyscale(spec)[chs], y=ytransform.(specdata[i][chs]), Geom.step, Theme(default_color=clr)))
    end
	maxE = ismissing(xmax) ? maxE : xmax
	append!(klms, autoklms ? mapreduce(s->NeXLSpectrum.elements(s, true, []), append!, specs) : [])
	if length(klms)>0
		tr(elm::Element) = characteristic(elm, alltransitions, 1.0e-3, maxE)
	    tr(cxr::CharXRay) = [ cxr ]
	    pklms = mapreduce(klm->tr(klm),append!,klms)
		if length(pklms)>0
			push!(layers, klmLayer(specdata,pklms))
		end
	end
	if length(edges)>0
		shs(elm::Element) = atomicshells(elm, maxE)
		shs(ash::AtomicShell) = [ ash ]
		pedges = mapreduce(ash->shs(ash), append!, edges)
		if length(pedges)>0
			push!(layers, edgeLayer(0.5*maxI,pedges))
		end
	end
    plot(layers...,
        Guide.XLabel("Energy (eV)"), Guide.YLabel(ylbl),
        Scale.x_continuous(format=:plain), Scale.y_continuous(format=:plain),
		Guide.manual_color_key(length(specs)>1 ? "Spectra" : "Spectrum", names, colors),
        Coord.Cartesian(ymin=0, ymax=ytransform(yscale*maxI), xmin=xmin, xmax=maxE))
end

"""
    plot(fd::FilteredData)

Plot the fitting filter and the spectrum from which it was derived.  Vertical red
lines represent the principle ROC.
"""
function Gadfly.plot(fd::FilteredDatum)
    p1=layer(x=fd.ffroi, y=fd.filtered, Geom.step, Theme(default_color=NeXLPalette[1]),
        xintercept=[fd.roi.start, fd.roi.stop], Geom.vline(color=["firebrick","firebrick"]))
    p2=layer(x=fd.ffroi, y=fd.data, Geom.step, Theme(default_color=NeXLPalette[2]))
	lyrs = ( p1, p2 )
    lbls = [ "Filtered", "Raw" ]
	if (fd isa FilteredReference) && (!ismissing(fd.back))
        bck = fd.data[fd.roi.start-fd.ffroi.start+1:fd.roi.stop-fd.ffroi.start+1]-fd.back
		p3=layer(x=fd.roi, y=bck, Geom.step, Theme(default_color=NeXLPalette[3]))
		lyrs = ( p1, p2, p3 )
        push!(lbls, "Background")
	end
    plot(lyrs...,Coord.cartesian(xmin=fd.ffroi.start, xmax=fd.ffroi.stop),
            Guide.XLabel("Channels"), Guide.YLabel("Counts"), Guide.title(repr(fd)),
            Guide.manual_color_key("", lbls, NeXLPalette[1:length(lbls)]))
end

"""
    Gadfly.plot(ffr::FilterFitResult)

Plot the unknown and the residual spectrum.
"""
function Gadfly.plot(ffr::FilterFitResult, roi::Union{Missing, UnitRange{Int}}=missing)
	if ismissing(roi)
		roi = ffr.roi
	end
	p1=layer(x=roi, y=ffr.raw[roi], Geom.step, Theme(default_color=NeXLPalette[1]))
    p2=layer(x=roi, y=ffr.residual[roi], Geom.step, Theme(default_color=NeXLPalette[2]))
    plot(p1,p2,Coord.cartesian(xmin=roi.start, xmax=roi.stop),
            Guide.XLabel("Channels"), Guide.YLabel("Counts"), Guide.title("$(ffr.label)"))
end
