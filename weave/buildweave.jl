using Weave

let start_dir = pwd()
    cd(@__DIR__)
    outpath = normpath(joinpath(@__DIR__, "..", "docs", "build"))
    @show outpath
    if !isdirpath(outpath)
        mkpath(outpath)
    end

    weave("errorbars.jmd", out_path=joinpath(outpath,"errorbars.html"))
    weave("K412fit.jmd", out_path=joinpath(outpath,"K412fit.html"))
    # weave("quantAMDglass.jmd", out_path=joinpath(outpath,"quantAMDglass.html"))
    weave("XRFspectra.jmd", out_path=joinpath(outpath,"XRFSpectra.html"))
    weave("K412quick.jmd", out_path=joinpath(outpath,"K412quick.html"))

    cd(start_dir)
end
