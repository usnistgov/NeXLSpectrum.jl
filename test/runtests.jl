using Test
using NeXLSpectrum

using DataDeps

register(DataDep("map[15]Artifact",
    """
    Dataset:    HyperSpectrum Test Data Set
    Author:     Nicholas W. M. Ritchie
    License:    CC-SA 3.0
    """,
    "https://drive.google.com/uc?export=download&id=1C93kn9-EIXXMDPcqJ9E4Xt4j9qfs5eeX",
    "c38ea6be06c6ef7eba427c50ecd44bc7dc4c63b5f8ef21897a8a51771edb7473",
    post_fetch_method=DataDeps.unpack
))
ENV["DATADEPS_ALWAYS_ACCEPT"]="true"

include("window.jl")
include("detefficiency.jl")
include("spectrum.jl")
include("emsa.jl")
include("brukerspx.jl")
include("filterfit.jl")
include("hyperspectrum.jl")
include("line.jl")
include("multikev.jl")
include("standardize.jl")
include("multidet.jl")
