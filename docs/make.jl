using Documenter, NeXLSpectrum, Gadfly, Colors

weavedocs = ( "continuummodel", "errorbars", "K412fit", "K412quick", "XRFSpectra" )

rebuildweave = true # !all(map(weavedocs) do name
#    isfile(joinpath(@__DIR__, "src", "$name.md"))
#end)

if rebuildweave
    map(name->rm(joinpath(@__DIR__, "src","$name.md")), filter(wd->isfile(joinpath(@__DIR__, "src","$name.md")),weavedocs))
    let curr=pwd()
        try
            include(joinpath(@__DIR__, "..","weave","buildweave.jl"))
        finally
            cd(curr)
        end
    end
end

pages = [
    "Home" => "index.md",
    "Spectrum Methods" => "spectrum.md",
    "HyperSpectrum Methods" => "hyperspectrum.md",
    "Fitting K412 (simple API)" => "k412refs.md",
    "Fitting K412 (flexible API)" => "K412fit.md",
    "Fitting K412 (quick fit)" => "K412quick.md",
    # "Quantifying AMD-6005a glass" => "quantAMDglass.md",
    "Fitting XRF Spectra" => "XRFspectra.md",
    "Lovely Error Bars" => "errorbars.md",
    "Modeling the Continuum" => "continuummodel.md",
    "Methods" => "methods.md",
    "Precompilation" => "precompile.md"
 ]

makedocs(
    modules = [NeXLSpectrum],
    sitename = "NeXLSpectrum.jl",
    pages = pages,
    checkdocs=:exports
)

function addNISTHeaders(htmlfile::String)
    # read HTML
    html = transcode(String,read(htmlfile))
    # Find </head>
    i = findfirst(r"</[Hh][Ee][Aa][Dd]>", html)
    # Already added???
    j = findfirst("nist-header-footer", html)
    if isnothing(j) && (!isnothing(i))
        # Insert nist-pages links right before </head>
        res = html[1:i.start-1]*
            "<link rel=\"stylesheet\" href=\"https://pages.nist.gov/nist-header-footer/css/nist-combined.css\">\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/jquery-1.9.0.min.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/nist-header-footer.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            html[i.start:end]
        write(htmlfile, res)
        println("Inserting NIST header/footer into $htmlfile")
    end
    return htmlfile
end

addNISTHeaders(joinpath(@__DIR__, "build","index.html"))
addNISTHeaders.(map(name->joinpath(@__DIR__, "build", splitext(name)[1], "index.html"), map(p->p.second, pages[2:end])))

#map(name->rm(joinpath("src","$(splitext(name)[1]).md")), weavedocs)

# deploydocs(repo = "github.com/NeXLSpectrum.jl.git")
