using Compose, Colors, CSV
using DataFrames


"""
    add_acquisition_points(imgpath::AbstractString, mime::AbstractString, x::Real, y::Real, fovx::Real, fovy::Real, pts::DataFrame)

A generic function for adding marker points to an image.  The image at `imgpath` of MIME-type `mime` centered at stage position `x`,`y`
with X and Y field-of-views `fovx`, `fovy` is annotated.  The point markers are placed at the coordinates pts[:,:X], pts[:,:Y].
"""
function add_acquisition_points(
    img::AbstractString,
    mime::AbstractString,
    x::Real,
    y::Real,
    fovx::Real,
    fovy::Real,
    pts::DataFrame,
)
    col = RGBA(colorant"yellow", 0.4) # Transparent yellow
    x0, y0 = x + 0.5 * fovx, y - 0.5 * fovy
    xx, yy = ((x0 .- pts[:, :X]) ./ fovx), (pts[:, :Y] .- y0) ./ fovy
    onimg =
        filter(i -> xx[i] > 0.0 && xx[i] < 1.0 && yy[i] > 0.0 && yy[i] < 1.0, eachindex(xx))
    compose(
        context(units = UnitBox(0.0, 0.0, 1.0, 1.0)),
        (
            context(),
            circle(xx[onimg], yy[onimg], [0.005]),
            Compose.stroke(col),
            linewidth(0.003w),
            fill("transparent"),
        ),
        (context(), bitmap(mime, read(img), 0.0, 0.0, 1.0, 1.0)),
    )
end

"""
    add_annotations(imgpath::AbstractString, modes::Symbol...; color = RGBA(colorant"yellow", 0.4), mime::AbstractString = "image/png")

Adds annotations to a SEMAntics image from the `images.txt` and `spectra.txt` meta-data files.

  Modes
  * `:ScaleBar`           Add a scale bar
  * `:AcquisitionPoints`  Mark all spectrum acquisition points.
"""
function add_annotations(
    imgpath::AbstractString,
    modes::Symbol...;
    color = RGBA(colorant"yellow", 0.4),
    mime::AbstractString = "image/png",
)
    imgdir, imgfn = splitdir(imgpath)
    df = CSV.read(joinpath(imgdir, "images.txt")) |> DataFrame
    img = splitext(imgfn)[1][1:end-3]
    rws = filter(r -> r[:Name] == img, df)
    if nrow(rws) > 0
        r = last(rws)
        aspect, fov = r[:YDim] / r[:XDim], r[:FOV]
        contexts = Any[context(units = UnitBox(0.0, 0.0, 1.0, aspect))]
        if :ScaleBar in modes
            sc = 10.0^floor(log10(fov))
            scale =
                fov / sc > 8.0 ? 5.0 * sc :
                (fov / sc > 4.0 ? 2.0 * sc : (fov / sc > 2.0 ? 1.0 * sc : 0.5 * sc))
            backcol = RGBA(0.8, 0.8, 0.8, 0.1)
            unit = (sc < 0.001 ? "nm" : (sc < 1.0 ? "μm" : "mm"))
            val = sc < 0.001 ? scale * 1.0e6 : (sc < 1.0 ? scale * 1.0e3 : scale)
            push!( # Scale bar
                contexts,
                (
                    context(),
                    rectangle(0.96 - scale / fov, 0.036, scale / fov, 0.006),
                    fill(color),
                ),
            )
            push!( # Text for scale bar
                contexts,
                (
                    context(),
                    text(0.96 - 0.5 * scale / fov, 0.06, "$val $unit", hcenter, vtop),
                    fontsize(0.04h),
                    font("Helvetica"),
                    fill(color),
                ),
            )
            push!( # Transparent gray background
                contexts,
                (
                    context(),
                    rectangle(0.954 - scale / fov, 0.03, scale / fov + 0.012, 0.07),
                    fill(backcol),
                ),
            )
        end
        if :AcquisitionPoints in modes
            x0, y0 = r[:X] + 0.5 * fov, r[:Y] - 0.5 * fov * aspect
            fovx, fovy = r[:FOV], aspect * r[:FOV]
            pts = CSV.read(joinpath(imgdir, "spectra.txt")) |> DataFrame
            xx, yy = ((x0 .- pts[:, :X]) ./ fovx), (pts[:, :Y] .- y0) ./ fovy
            onimg = filter(
                i -> xx[i] > 0.0 && xx[i] < 1.0 && yy[i] > 0.0 && yy[i] < 1.0,
                eachindex(xx),
            )
            push!( # X markers
                contexts,
                (
                    context(),
                    xgon(xx[onimg], yy[onimg], [0.02 * w], [4]),
                    stroke(color),
                    linewidth(0.003 * w),
                    fill("transparent"),
                ),
            )
        end
        # Image
        push!(contexts, (context(), bitmap(mime, read(imgpath), 0.0, 0.0, 1.0, aspect)))
        return compose(contexts...)
    end
end

if false
    path = "C:\\Users\\nritchie\\Desktop\\Topsham Maine Uraninite Hex Crystal"

    df = vcat(
        (CSV.read(joinpath(path, "Edge $i", "spectra.txt")) |> DataFrame for i = 1:6)...,
        (CSV.read(joinpath(path, "Traverse $i", "spectra.txt")) |> DataFrame for i = 1:3)...,
    )

    set_default_graphic_size(20cm, 20cm)

    add_acquisition_points(
        joinpath(path, "Manual", "Macro[0].png"),
        "image/png",
        -0.6720,
        -0.02550,
        7.0,
        7.0,
        df,
    ) |> SVG(joinpath(path, "annot macro.svg"))
    add_acquisition_points(
        joinpath(path, "Manual", "Macro2[0].png"),
        "image/png",
        -0.5700,
        -0.6950,
        14.0,
        14.0,
        df,
    ) |> SVG(joinpath(path, "annot macro2.svg"))
    add_annotations(joinpath(path, "Edge 1", "Point 10[0].png"), :ScaleBar, :AcquisitionPoints)
end