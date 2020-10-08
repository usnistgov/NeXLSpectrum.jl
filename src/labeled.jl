using Compose
using FileIO
using ImageIO

struct LabeledView
    label::String
    image::Array
end


"""
    labeledimage(label::String, image::Array)

Displays a string label below an image.  Useful in Weave, Juno or Jupyter.
"""
function labeledimage(label::String, image::Array)
    io = IOBuffer(maxsize=10*1024*1024)
    save(Stream(format"PNG",io), image)
    pix = max(size(image,1),size(image,2))
    scaleX, scaleY = size(image,1)/pix, size(image,2)/pix
    return compose(context(),
            (context(0.5*(1.0-0.8*scaleX), 0.5*(1.0-0.8*scaleY), 0.8*scaleX, 0.8*scaleY), bitmap("image/png", take!(io), 0.0, 0.0, scaleX, scaleY)),
            (context(0.1, 0.8, 0.8, 0.2), text(0.5, 0.9, label, hcenter, vbottom), stroke("transparent"), fill("black"), fontsize(0.06h))
    )
end



"""
    labeledimages(labels::AbstractVector{<:AbstractString}, images::AbstractVector{<:AbstractArray}; ncols=3, halign = hleft)

Create a matrix of labeled images.  Useful in Weave, Juno or Jupyter.
"""
function labeledimages(labels::AbstractVector{<:AbstractString}, images::AbstractVector{<:AbstractArray}; ncols=3, halign = hleft)
    @assert length(labels)==length(images)
    nrows = (length(images)+ncols-1) ÷ ncols
    sc = 1.0/max(ncols, nrows)
    tmp = []
    for i in eachindex(labels)
        r, c = (i-1) ÷ ncols, (i-1) % ncols
        push!(tmp, (context(c*sc, r*sc, sc, sc), labeledimage(labels[i], images[i])))
    end
    return compose(context(), tmp...)
end
