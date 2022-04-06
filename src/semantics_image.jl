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
        # Single class to replace fails in 1.6
        rg = Regex(reduce((a,b)->replace(a,b), ( "["=>"\\[", "]"=>"\\]", "_" => "[.\\\\:_]"); init=string(nm[1:end-3])))
        mimgs = filter(r -> !isnothing(match(rg, r[:Name])), imgs)
        if nrow(mimgs) > 0
            # Since an image can be written more than once, take the last...
            x, y, fov = mimgs[end, :X], mimgs[end, :Y], mimgs[end, :FOV]
            pix = fov / size(img, 2) # square pixels determined by horizontal FOV
            img = AxisArray(
                img, #
                Axis{:y}(StepRangeLen(y, pix, size(img,1), size(img,1)÷2)),
                Axis{:x}(StepRangeLen(x, -pix, size(img,2), size(img,2)÷2))
            )
                    
        else
            @info "There is no meta-data in images.txt for images matching $rg for $fn."
        end
    else
        @info "There is no meta-data file images.txt for images matching $rg for $fn."
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

function build_scale_bar(img::AxisArray; marker = HSVA(60, 1, 1, 0.3))
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

function build_acquisition_points(
    img::AxisArray,
    coords::AbstractArray{<:AbstractDict{Symbol, <:AbstractFloat}};
    label = true, 
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
        push!(items, (context(), Compose.text(fx.(x), y, repr.(eachindex(coords)), fh.(x), [ vcenter]), fill(marker), stroke("transparent")))
    end
    return items
end

function build_thumbnails(
    img::AxisArray, 
    thumbnails::Union{Nothing, AbstractArray{<:AxisArray}}=nothing;
    label = true, 
    marker = HSVA(60, 1, 1, 0.3)
)
    ya, xa = axisvalues(img)[1:2]
    rx, ry, rw, rh = Float64[], Float64[], Float64[], Float64[]
    for thumb in thumbnails
        yt, xt = axisvalues(thumb)[1:2]
        y1o = (yt[1]-ya[1])/(ya[end]-ya[1])
        yeo = (yt[end]-ya[1])/(ya[end]-ya[1])
        x1o = (xt[1]-xa[1])/(xa[end]-xa[1])
        xeo = (xt[end]-xa[1])/(xa[end]-xa[1])
        if (max(x1o,xeo)>0.0) && (min(x1o,xeo)<1.0) &&
            (max(y1o,yeo)>0.0) && (min(y1o,yeo)<1.0)
            push!(rx, x1o)
            push!(ry, y1o)
            push!(rw, xeo-x1o)
            push!(rh, yeo-y1o)
        end
    end
    #  Place rectangles at image points
    items = Tuple[ (context(), Compose.rectangle(rx, ry, rw, rh), fill(RGBA(0.1, 0.1, 0.1, 0.3)), stroke(marker)) ]
    if label
        push!(items, (context(), Compose.text(rx .+ 0.001, ry .+ 0.005, repr.(eachindex(thumbnails)), [ hleft ], [ vtop ]), fill(marker), stroke("transparent")))
    end
    return items
end

"""
    annotate(
        img::AxisArray; 
        scale::Bool = true,
        coords::Union{Nothing, AbstractArray{<:AbstractDict{Symbol, <:AbstractFloat}}} = nothing,
        spectra::Union{Nothing, <:AbstractArray{<:Spectrum}} = nothing,
        thumbnails::Union{Nothing, AbstractArray{<:AxisArray}}=nothing,
        labelcoords::Bool = true,
        labelthumbnails::Bool = true,
        marker::Colorant = HSVA(60, 1, 1, 0.3) 
    )

Add annotations like scale-bars, acquisition points, image areas to an image
and display the result.

  * `coords` and `spectra` display as circles (with/without numeric labels).
  * `thumbnails` display as rectangles overlaying the image
  * `scale` displays as a scale-bar in the upper-right corner
  * `marker` is the color of the annotations which is by-default slightly transparent
"""
function annotate(
    img::AxisArray; 
    scale::Bool = true,
    coords::Union{Nothing, AbstractArray{<:AbstractDict{Symbol, <:AbstractFloat}}} = nothing,
    spectra::Union{Nothing, <:AbstractArray{<:Spectrum}} = nothing,
    thumbnails::Union{Nothing, AbstractArray{<:AxisArray}}=nothing,
    labelcoords::Bool = true,
    labelthumbnails::Bool = true,
    marker::Colorant = HSVA(60, 1, 1, 0.3) 
)
    items = Tuple[]
    if scale
        append!(items, build_scale_bar(img, marker=marker))
    end
    allcoords = AbstractDict{Symbol, <:AbstractFloat}[]
    if !isnothing(coords)
        append!(allcoords, coords)
    end
    if !isnothing(spectra)
        foreach(sp->push!(allcoords, sp[:StagePosition]), spectra)
    end
    if length(allcoords)>0
        append!(items, build_acquisition_points(img, allcoords; label=labelcoords, marker=marker))
    end
    if !isnothing(thumbnails)
        append!(items, build_thumbnails(img, thumbnails; label=labelthumbnails, marker=marker))
    end
    io = IOBuffer()
    ImageIO.save(Stream{format"PNG"}(io), img)
    return compose(
        context(0.0, 0.0, 1.0, 1.0),
        items...,
        ( context(), Compose.bitmap("image/png", take!(io), 0.0, 0.0, 1.0, 1.0) ),
    )
end

"""
    shannon_entropy(img::AbstractArray{Gray{N0f8}})

Computes the log2-entropy of the data in the image.
The entropy(...) in Images.jl 24.1 is buggy and is removed in 25.0
"""
function shannon_entropy(img::AbstractArray{Gray{N0f8}})
    # c = @MVector zeros(Int(typemax(UInt8))+1) # Reduces allocations to zero..
    c = zeros(Int(typemax(UInt8))+1) 
    for v in img
        @inbounds c[reinterpret(UInt8, v.val)+1] += 1.0/length(img)
    end
    l2(v) = v==0 ? 0.0 : log2(v)
    return -sum(v->v*l2(v), c)
end