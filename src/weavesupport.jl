import Compose
import .Weave
try
    import Cairo
catch
    @warn("Cairo.jl is required to be installed to generate raster images")
end

Base.showable(m::MIME"application/svg", ctx::Context) = true
Base.showable(m::MIME"application/png", ctx::Context) = true

function Base.display(report::Weave.Report, m::MIME"image/svg+xml", ctx::Context)
    chunk = report.cur_chunk

    w = chunk.options[:fig_width]Compose.inch
    h = chunk.options[:fig_height]Compose.inch
    format = chunk.options[:fig_ext]
    dpi = chunk.options[:dpi]

    full_name, rel_name = Weave.get_figname(report, chunk, ext = format)

    push!(report.figures, rel_name)
    report.fignum += 1
    if format == ".svg"
        Compose.draw(Compose.SVG(full_name, w, h), ctx)
    elseif format == ".png"
        Compose.draw(Compose.PNG(full_name, w, h, dpi = dpi), ctx)
    elseif format == ".pdf"
        Compose.draw(Compose.PDF(full_name, w, h), ctx)
    elseif format == ".ps"
        Compose.draw(Compose.PS(full_name, w, h), ctx)
    elseif format == ".tex"
        Compose.draw(Compose.PGF(full_name, w, h, true), ctx)
    else
        @warn("Can't save figure. Unsupported format, $format")
    end
end

@info "Weave support loaded into NeXLSpectrum"
