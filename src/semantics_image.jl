using FileIO, CSV, ImageAxes, Compose

"""
    readSEManticsImage(fn::AbstractString)

Read a SEMantics PNG image and the create an `ImageAxes` (`AxisArray`) which associates
`:x`, `:y` coordinates with the image axes according to the X, Y and FOV in the
images.txt file in the same directory.
"""
function readSEManticsImage(fn::AbstractString)
    img, sp = load(fn), splitpath(fn)
    img_txt = joinpath(sp[1:end-1]..., "images.txt")
    if isfile(img_txt)
        imgs = CSV.File(img_txt, delim = "\t") |> DataFrame
        (nm, _) = splitext(sp[end])
        rg = Regex(replace(nm[1:end-3], "_" => "[.\\\\:]"))
        mimgs = filter(r -> !isnothing(match(rg, r[:Name])), imgs)
        if nrow(mimgs) > 0
            # Since an image can be written more than once, take the last...
            x, y, fov = mimgs[end, :X], mimgs[end, :Y], mimgs[end, :FOV]
            pix = fov / size(img, 2) #square pixels
            img = AxisArray(
                img, #
                Axis{:y}(StepRangeLen(y, pix, size(img,1), size(img,1)÷2)),
                Axis{:x}(StepRangeLen(x, -pix, size(img,2), size(img,2)÷2))
            )
                    
        else
            @info "There is no meta-data in image.txt for images matching $rg."
        end
    else
        @info "There is no meta-data file images.txt for images matching $rg."
    end
    return img
end

function _ns_closest(srl::AbstractRange, val)
    idx = searchsortedfirst(srl, val, rev = step(srl) < 0.0)
    ip = idx + (step(srl) < 0.0 ? -1 : 1)
    return (val < min(srl[1], srl[end])) || (val > max(srl[1], srl[end])) ? nothing : #
           (abs(srl[idx] - val) < abs(srl[ip] - val) ? idx : ip)
end

function _ns_closest(img::AxisArray, coord)
    return ( 
        _ns_closest(ImageAxes.axisvalues(img)[2], coord[1]), 
        _ns_closest(ImageAxes.axisvalues(img)[1], coord[2]), 
    )
end

function scale_bar(img::AxisArray; marker = HSVA(60, 1, 1, 0.3))
    ax = AxisArrays.axisvalues(img)[1]
    fov = abs(ax[end] - ax[1]) / 2.0
    sc = 10.0^floor(log10(fov))
    l = (fov / sc >= 5.0 ? 5.0 : (fov / sc >= 2.5 ? 2.5 : 1.0)) * sc
    dp = l / abs(ax[end] - ax[1])
    xo, yo = 0.02, 0.02
    return Tuple[
        (context(), Compose.line([(1.0 - 2.0 * xo - dp, yo + 0.07), (1.0 - 2.0 * xo, yo + 0.07)]), stroke(marker), linewidth(1.0)),
        (context(), Compose.text(1.0 - 2.0 * xo - dp / 2.0, yo + 0.045, "$(round(Int,1000.0*l)) μm", hcenter, vbottom), fill(marker), stroke("transparent"), fontsize(8.0)),
        (context(), Compose.rectangle(1.0 - 3.0 * xo - dp, yo, dp + 2.0 * xo, yo + 0.07), stroke("transparent"), fill(RGBA(0.1, 0.1, 0.1, 0.3)))
    ]
end

"""
    mark_acquisition_points(
        img::AxisArray, 
        coords::AbstractArray{<:AbstractDict{Symbol,<:AbstractFloat}}; 
        label = true, 
        scale = true, 
        marker = HSVA(60, 1, 1, 0.3)
    )
    mark_acquisition_points(
        img::AxisArray, 
        specs::AbstractArray{<:Spectrum}; 
        label = true, 
        scale = true, 
        marker = HSVA(60, 1, 1, 0.3)
    )

Mark the locations on the ImageAxes scaled image from which the list of Spectrum objects
was collected.  Optionally, place numeric index labels and a scale marker on the image.

`coords` contains an array of dictionaries with keys `:X` and `:Y` with the x- and y-
coordinates.

This function assumes that the spectra were collected from the stage coordinates in the
`sp=spec[:StagePosition]` property as `sp[:X]` and `sp[:Y]`. 
"""
function mark_acquisition_points(
    img::AxisArray, 
    coords::AbstractArray{<:AbstractDict{Symbol, <:AbstractFloat}}; 
    label = true, 
    scale = true, 
    marker = HSVA(60, 1, 1, 0.3)
)
    x, y = Float64[], Float64[]
    for coord in coords
        pt = _ns_closest(img, (coord[:X], coord[:Y]))
        if !any(isnothing.(pt))
            push!(x, pt[1] / size(img, 1))
            push!(y, pt[2] / size(img, 2))
        end
    end
    #  Place circles at acquisition point
    items = Tuple[(context(), Compose.circle(x, y, [0.005]), fill(marker), stroke("transparent"))]
    if label
        fx(x) = x < 0.9 ? x + 0.01 : x - 0.01
        fh(x) = x < 0.9 ? hleft : hright
        push!(items, (context(), Compose.text(fx.(x), y, repr.(eachindex(coords)), fh.(x), [vcenter]), fill(marker), stroke("transparent")))
    end
    # Scale bar
    if scale
        append!(items, scale_bar(img, marker = marker))
    end
    # Draw the image to a context()
    io = IOBuffer()
    ImageIO.save(Stream{format"PNG"}(io), img)
    return compose(
        context(0.0, 0.0, 1.0, 1.0),
        items...,
        (context(), Compose.bitmap("image/png", take!(io), 0.0, 0.0, 1.0, 1.0)),
    )
end

function mark_acquisition_points(
    img::AxisArray, 
    specs::AbstractArray{<:Spectrum}; 
    kwargs...
)
    coords = [ spec[:StagePosition] for spec in specs ]
    mark_acquisition_points(img, coords; kwargs...)
end
