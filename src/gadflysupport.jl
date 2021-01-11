using .Gadfly
using Colors
using Printf

const NeXLSpectrumStyle = style(
    background_color = nothing,
    panel_fill = RGB(253 / 255, 253 / 255, 241 / 255),
    grid_color = RGB(255 / 255, 223 / 255, 223 / 255),
    grid_color_focused = RGB(255 / 255, 200 / 255, 200 / 255),
    grid_line_style = :solid,
    major_label_color = RGB(32 / 255, 32 / 255, 32 / 255),
    major_label_font_size = 9pt,
    panel_stroke = RGB(32 / 255, 32 / 255, 32 / 255),
    plot_padding = [2pt],
    key_title_font_size = 9pt,
    key_position = :right, # :bottom
)


"""
    plot(
	    specs::Union{Spectrum...,AbstractVector{Spectrum{<:Real}}};
	    klms=[],
	    edges=[],
		escapes=[],
		coincidences=[],
	    autoklms = false,
	    xmin=0.0,
	    xmax=missing,
	    norm=:None,
	    yscale=1.05,
	    ytransform = identity,
		style=NeXLSpectrumStyle,
		palette=NeXLPalette
    )::Plot

Required:

	specs::AbstractVector{Spectrum};

Named:

    klms = [ Element &| CharXRay ]
	edges = [ Element &| AtomicSubShell ]
	escapes = [ CharXRay ],
	coincidences = [ CharXRay ]
	autoklms = false # Add KLMs based on elements in spectra
	xmin = 0.0 # Min energy (eV)
	xmax = missing # Max energy (eV) (defaults to max(:BeamEnergy))
	norm = NoScaling() | ScaleDoseWidth() | ScaleDose() | ScaleSum() | ScaleROISum() | ScalePeak() | (<: SpectrumScaling)()
	yscale = 1.05 # Fraction of max intensity for ymax over [max(lld,xmin):xmax]
	ytransform = identity | log10 | sqrt | ??? # How to transform the counts data before plotting
	style=NeXLSpectrumStyle (or another Gadfly.style)
	palette = NeXLPalette | Colorant[ ... ] # Colors for spectra...
    customlayers = Gadfly.Layer[] # Allows additional plot layers to be added

Plot a multiple spectra on a single plot using Gadfly.

    Gadfly.plot(
      ffr::FilterFitResult,
      roi::Union{Missing,UnitRange{Int}} = missing;
      palette = NeXLPalette,
      style = NeXLSpectrumStyle, xmax=missing)

Plot the spectrum, residual and k-ratios using Gadfly.

    Gadfly.plot(fr::FilteredReference)

Plot a filtered reference spectrum.
"""
Gadfly.plot( #
    specs::AbstractVector{Spectrum{<:Real}};
    klms = [],
    edges = [],
    escapes = [],
    coincidences = [],
    autoklms = false,
    xmin = 0.0,
    xmax = missing,
    norm = :None,
    yscale = 1.05,
    ytransform = identity,
    style = NeXLSpectrumStyle,
    palette = NeXLPalette,
)::Plot = plot( #
    specs...,
    klms = klms,
    edges = edges,
    escapes = escapes,
    coincidences = coincidences,
    autoklms = autoklms,
    xmin = xmin,
    xmax = xmax,
    norm = norm,
    yscale = yscale,
    ytransform = identity,
    style = style,
    palette = palette,
)

function Gadfly.plot(
    specs::Spectrum{<:Real}...;
    klms = [],
    edges = [],
    escapes = [],
    coincidences = [],
    autoklms = false,
    xmin = 0.0,
    xmax = missing,
    legend=true,
    norm = NoScaling(),
    yscale = 1.05,
    ytransform = identity,
    style = NeXLSpectrumStyle,
    palette = NeXLPalette,
    customlayers = Gadfly.Layer[]
)::Plot
    function klmLayer(specdata, cxrs::AbstractArray{CharXRay})
        d = Dict{Any,Vector{CharXRay}}()
        for cxr in cxrs
            d[(element(cxr), shell(cxr))] = push!(get(d, (element(cxr), shell(cxr)), []), cxr)
        end
        x, y, label = [], [], []
        for cs in values(d)
            br = brightest(cs)
            ich = maximum(specdata[i][channel(energy(br), spec)] for (i, spec) in enumerate(specs))
            if ich > 0
                for c in cs
                    push!(x, energy(c))
                    push!(y, ytransform(ich * weight(c)))
                    push!(label, weight(c) > 0.1 ? "$(element(c).symbol)" : "")
                end
            end
        end
        return layer(
            x = x,
            y = y,
            label = label,
            Geom.hair,
            Geom.point,
            Geom.label(position = :above),
            Theme(default_color = colorant"gray"),
        )
    end
    function edgeLayer(maxI, ashs::AbstractArray{AtomicSubShell})
        d = Dict{Any,Vector{AtomicSubShell}}()
        for ash in ashs
            d[(element(ash), shell(ash))] = push!(get(d, (element(ash), shell(ash)), []), ash)
        end
        x, y, label = [], [], []
        for ass in values(d)
            br = ass[findmax(capacity.(ass))[2]]
            for ash in ass
                push!(x, energy(ash))
                push!(y, ytransform(maxI * capacity(ash) / capacity(br)))
                push!(label, "$(ash)")
            end
        end
        return layer(
            x = x,
            y = y,
            label = label,
            Geom.hair,
            Geom.label(position = :right),
            Theme(default_color = colorant"lightgray"),
        )
    end
    function siEscapeLayer(crxs)
        x, y, label = [], [], []
        for xrs in crxs
            eesc = energy(xrs) - energy(n"Si K-L3")
            if eesc > 0.0
                ich = maximum(get(specdata[i], channel(eesc, spec), 0.0) for (i, spec) in enumerate(specs))
                push!(x, eesc)
                push!(y, ytransform(ich))
                push!(label, "$(element(xrs).symbol)\nesc")
            end
        end
        return layer(
            x = x,
            y = y,
            label = label,
            Geom.hair,
            Geom.label(position = :above),
            Theme(default_color = colorant"black"),
        )
    end
    function sumPeaks(cxrs)
        x, y, label = [], [], []
        for (i, xrs1) in enumerate(cxrs)
            for xrs2 in cxrs[i:end] # place each pair once...
                eesc = energy(xrs1) + energy(xrs2)
                ich = maximum(get(specdata[i], channel(eesc, spec), 0.0) for (i, spec) in enumerate(specs))
                if ich > 0.0
                    push!(x, eesc)
                    push!(y, ytransform(ich))
                    push!(label, "$(element(xrs1).symbol)\n+\n$(element(xrs2).symbol)")
                end
            end
        end
        return layer(
            x = x,
            y = y,
            label = label,
            Geom.hair,
            Geom.label(position = :above),
            Theme(default_color = colorant"gray"),
        )
    end
    @assert length(specs)<=length(palette) "The palette must specify at least as many colors as spectra."
    specdata = [ scaledcounts(norm,s) for s in specs]
    ylbl = repr(norm)
    maxI, maxE, maxE0 = 16, 1.0e3, 1.0e3
    names, layers = String[], Layer[]
    append!(layers, customlayers)
    for (i, spec) in enumerate(specs)
        mE = ismissing(xmax) ? get(spec, :BeamEnergy, energy(length(spec), spec)) : convert(Float64,xmax)
        mE0 = get(spec, :BeamEnergy, missing)
        chs = max(1, channel(convert(Float64,xmin), spec)):min(length(spec),channel(mE, spec))
        mchs = max(channel(200.0,spec), chs[1], lld(spec)):min(length(specdata[i]),chs[end])  # Ignore zero strobe...
        maxI = max(maxI, maximum(specdata[i][mchs]))
        maxE = max(maxE, mE)
        maxE0 = ismissing(mE0) ? maxE : max(maxE, mE0)
        push!(names, spec[:Name])
        ly = Gadfly.layer(x = energyscale(spec)[chs], y = ytransform.(specdata[i][chs]), #
            Geom.step, Theme(default_color = palette[i]))
        append!(layers, ly)
    end
    append!(klms, autoklms ? mapreduce(s -> elms(s, true, []), union!, specs) : [])
    if length(klms) > 0
        tr(elm::Element) = characteristic(elm, alltransitions, 1.0e-3, maxE0)
        tr(cxr::CharXRay) = [cxr]
        pklms = mapreduce(klm -> tr(klm), append!, klms)
        if length(pklms) > 0
            append!(layers, klmLayer(specdata, pklms))
        end
    end
    if length(edges) > 0
        shs(elm::Element) = atomicsubshells(elm, maxE0)
        shs(ash::AtomicSubShell) = [ash]
        pedges = mapreduce(ash -> shs(ash), append!, edges)
        if length(pedges) > 0
            append!(layers, edgeLayer(0.5 * maxI, pedges))
        end
    end
    if length(escapes) > 0
        append!(layers, siEscapeLayer(escapes))
    end
    if length(coincidences) > 0
        append!(layers, sumPeaks(coincidences))
    end
    Gadfly.with_theme(style) do
        leg = legend ? 
            tuple( Guide.manual_color_key(length(specs) > 1 ? "Spectra" : "Spectrum", names, palette[1:length(specs)]) ) : 
            tuple(  )
        try
            plot(
                layers...,
                Guide.XLabel("Energy (eV)"),
                Guide.YLabel(ylbl),
                Scale.x_continuous(format = :plain),
                Scale.y_continuous(format = :plain),
                Coord.Cartesian(ymin = 0, ymax = ytransform(yscale * maxI), xmin = convert(Float64,xmin), xmax = maxE),
                leg...,
            ) 
        catch
            plot(
                layers...,
                Guide.XLabel("Energy (eV)"),
                Guide.YLabel(ylbl),
                Scale.x_continuous(format = :plain),
                Scale.y_continuous(format = :plain),
                Coord.Cartesian(ymin = 0, ymax = ytransform(yscale * maxI), xmin = convert(Float64,xmin), xmax = maxE),
                leg...,
            )
        end
    end
end

function Gadfly.plot(ffr::FilterFitResult, roi::Union{Missing,UnitRange{Int}} = missing; palette = NeXLPalette, style = NeXLSpectrumStyle, xmax=missing)
    function defroi(ffrr) # Compute a reasonable default display ROI
        tmp = minimum( lbl.roi[1] for lbl in keys(ffrr.kratios)):maximum( lbl.roi[end] for lbl in keys(ffrr.kratios))
        return max(lld(ffr.label.spectrum),tmp[1]-length(ffrr.roi)÷40):min(tmp[end]+length(ffrr.roi)÷10,ffrr.roi[end])
    end
    roilt(l1, l2) = isless(l1.roi[1], l2.roi[1])
    roi = ismissing(roi) ? defroi(ffr) : roi
    layers = [
        layer(x = roi, y = ffr.residual[roi], Geom.step, Theme(default_color = palette[2])),
        layer(x = roi, y = ffr.raw[roi], Geom.step, Theme(default_color = palette[1])),
    ]
    miny, maxy, prev, i = minimum(ffr.residual[roi]), 3.0 * maximum(ffr.residual[roi]), -1000, -1
    for lbl in sort(collect(keys(ffr.kratios)), lt = roilt)
        if NeXLUncertainties.value(lbl, ffr) > 0.0
            # This logic keeps the labels on different lines (mostly...)
            i, prev = (lbl.roi[1] > prev + length(roi) ÷ 10) || (i == 6) ? (0, lbl.roi[end]) : (i + 1, prev)
            labels = ["", name(lbl.xrays)]
            # Plot the ROI
            push!(
                layers,
                layer(
                    x = [lbl.roi[1], lbl.roi[end]],
                    y = maxy * [0.4 + 0.1 * i, 0.4 + 0.1 * i],
                    label = labels,
                    Geom.line,
                    Geom.point,
                    Geom.label(position = :right),
                    Theme(default_color = "gray"),
                ),
            )
            # Plot the k-ratio as a label above ROI
            push!(
                layers,
                layer(
                    x = [0.5 * (lbl.roi[1] + lbl.roi[end])],
                    y = maxy * [0.4 + 0.1 * i],
                    label = [@sprintf("%1.4f", NeXLUncertainties.value(lbl, ffr))],
                    Geom.label(position = :above),
                    Theme(default_color = "gray"),
                ),
            )
        end
    end
    Gadfly.with_theme(style) do
        plot(
            layers...,
            Coord.cartesian(xmin = roi[1], xmax = ismissing(xmax) ? roi[end] : xmax, ymin = min(1.1 * miny, 0.0), ymax = maxy),
            Guide.XLabel("Channels"),
            Guide.YLabel("Counts"),
            Guide.title("$(ffr.label)"),
        )
    end
end

function Gadfly.plot(fr::FilteredReference; palette=NeXLPalette)
    roicolors = Colorant[ RGB(0.9, 1.0, 0.9), RGB(0.95, 0.95, 1.0)]
    layers=[
        layer(x=fr.ffroi, y=fr.data, Theme(default_color=palette[1]), Geom.step),
        layer(x=fr.ffroi, y=fr.filtered, Theme(default_color=palette[2]), Geom.step),
        layer(x=fr.roi, y=fr.charonly, Theme(default_color=palette[3]), Geom.step),
        layer(xmin=[fr.ffroi[1], fr.roi[1]], xmax=[fr.ffroi[end], fr.roi[end]], Geom.vband, color = roicolors ),
    ]
    try
        plot(layers..., Coord.cartesian(xmin=fr.ffroi[1]-length(fr.ffroi)÷10, xmax=fr.ffroi[end]+length(fr.ffroi)÷10),
            Guide.xlabel("Channel"), Guide.ylabel("Counts"), Guide.title(repr(fr.label)),
            Guide.manual_color_key("Legend",["Spectrum", "Filtered", "Char. Only", "Filter ROC", "Base ROC"], [ palette[1:3]..., roicolors...] ))
    catch
        plot(layers..., Coord.cartesian(xmin=fr.ffroi[1]-length(fr.ffroi)÷10, xmax=fr.ffroi[end]+length(fr.ffroi)÷10),
            Guide.xlabel("Channel"), Guide.ylabel("Counts"), Guide.title(repr(fr.label)),
        #    Guide.manual_color_key("Legend",["Spectrum", "Filtered", "Char. Only", "Filter ROC", "Base ROC"], [ palette[1:3]..., roicolors...] )
        )
        end
end

Gadfly.plot(ff::TopHatFilter, fr::FilteredReference) =
    spy(filterdata(ff, fr.ffroi), Guide.title(repr(fr.label)), Guide.xlabel("Channels"), Guide.ylabel("Channels"))

function Gadfly.plot(vq::VectorQuant, chs::UnitRange)
    colors = distinguishable_colors(
        size(vq.vectors,1)+2,
        [RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
        transform = deuteranopic,
    )[3:end]
    lyrs = mapreduce(i->layer(x=chs, y=vq.vectors[i,chs], Theme(default_color=colors[i]), Geom.line),append!,eachindex(vq.references))
    try
        plot(lyrs...,
            Guide.xlabel("Channel"), Guide.ylabel("Filtered"),
            Guide.manual_color_key("Vector", [ repr(r[1]) for r in vq.references ], color=Colorant[colors...]))
    catch
        plot(lyrs...,
            Guide.xlabel("Channel"), Guide.ylabel("Filtered"),
        #    Guide.manual_color_key("Vector", [ repr(r[1]) for r in vq.references ], color=Colorant[colors...])
        )
    end
end

function Gadfly.plot(deteff::DetectorEfficiency, emax=20.0e3)
    eff(ee) = efficiency(deteff, ee, π/2)
    plot(eff, 100.0, emax)
end

function plotandimage(plot::Gadfly.Plot, image::Array)
    io = IOBuffer(maxsize=10*1024*1024)
    save(Stream(format"PNG",io), image)
    pix = max(size(image,1),size(image,2))
    scaleX, scaleY = size(image,1)/pix, size(image,2)/pix
    return compose(context(),
        (context(0.0, 0.0, 0.8, 1.0), Gadfly.render(plot) ),
        (context(0.8, 0.0, 0.2, 1.0), bitmap("image/png", take!(io), 0.5*(1.0 - scaleX), 0.5*(1.0 - scaleY), scaleX, scaleY)),
    )
end

Gadfly.plot(ffp::FilterFitPacket; kwargs...) = plot( unique(spectra(ffp))...; klms=collect(elms(ffp)), kwargs...)

@info "Loading Gadfly support into NeXLSpectrum."
