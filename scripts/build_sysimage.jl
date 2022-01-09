# This script will build a `sysimage` for Windows.  
# Please submit a pull-request if you modify it to support additional operating systems.

# Load the libraries
using NeXLSpectrum, Gadfly, DataFrames 

# Load the package compiler
using PackageCompiler

# This specifies which libraries to include and also specifies a "precompile.jl" file which forces 
# the just-in-time compiler to compile frequently used functions.
PackageCompiler.create_sysimage(
    [ "Gadfly", "DataFrames", "BoteSalvatICX", "FFAST", "NeXLUncertainties", "NeXLCore", "NeXLMatrixCorrection", "NeXLSpectrum" ]; 
    sysimage_path=joinpath(homedir(), ".julia", "NeXLSysimage.dll"),
    precompile_execution_file=joinpath(pkgdir(NeXLSpectrum),"scripts", "precompile.jl"))
