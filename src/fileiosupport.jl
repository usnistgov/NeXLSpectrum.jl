# This file defines a mechanism for auto-sniffing, loading and saving spectrum files

abstract type SpectrumFileType end

function savespectrum(ty::Type{<:SpectrumFileType}, fio, data)
    @error "Saving to is not implemented for $ty. Probably never will be."
end

function loadspectrum(ty::Type{<:SpectrumFileType}, filename::AbstractString; kwargs...)
    badname(sp) =
        (!haskey(sp, :Name)) ||
        startswith(sp[:Name], "Bruker") ||
        startswith(sp[:Name], "Spectrum[")
    @assert isfile(filename) "$filename is not a file in loadspectrum(...)."
    return open(filename, read = true) do ios
        spec = loadspectrum(ty, ios; kwargs...)
        spec[:Filename] = filename
        spec[:Name] = badname(spec) ? splitext(splitpath(filename)[end])[1] : spec[:Name]
        return spec
    end
end

loadspectrum(
    ty::Type{<:SpectrumFileType},
    filename::AbstractString,
    det::EDSDetector;
    kwargs...,
) = apply(loadspectrum(ty, filename; kwargs...), det)

savespectrum(ty::Type{<:SpectrumFileType}, filename::AbstractString, spec::Spectrum) =
    open(filename, write = true) do ios
        savespectrum(ty, ios, spec)
        spec[:Filename] = filename
    end

# ---------------------------------------------------------------------------------------------------#

struct BrukerSPX <: SpectrumFileType end

loadspectrum(::Type{BrukerSPX}, ios::IO) = readbrukerspx(ios)
extensions(::Type{BrukerSPX}) = (".spx",)


sniff(::Type{BrukerSPX}, ios::IO) = detectbrukerspx(ios)

# ---------------------------------------------------------------------------------------------------#

struct BrukerPDZ <: SpectrumFileType end

loadspectrum(::Type{BrukerPDZ}, ios::IO) = readbrukerpdz(ios)
extensions(::Type{BrukerPDZ}) = (".pdz",)

sniff(::Type{BrukerPDZ}, ios::IO) = detectbrukerpdz(ios)

# ---------------------------------------------------------------------------------------------------#

struct ISOEMSA <: SpectrumFileType end

loadspectrum(::Type{ISOEMSA}, ios::IO, ty::Type{<:Real} = Float64) = readEMSA(ios, ty)
savespectrum(::Type{ISOEMSA}, ios::IO, spec::Spectrum) = writeEMSA(ios, spec)

extensions(::Type{ISOEMSA}) = (".msa", ".emsa", ".ems", ".txt")

sniff(::Type{ISOEMSA}, ios::IO) = isemsa(ios)

# ---------------------------------------------------------------------------------------------------#

struct ASPEXTIFF <: SpectrumFileType end

loadspectrum(::Type{ASPEXTIFF}, ios::IO; withImgs = false, astype::Type{<:Real} = Float64) =
    readAspexTIFF(ios, withImgs = withImgs, astype = astype)
extensions(::Type{ASPEXTIFF}) = (".tif",)

sniff(::Type{ASPEXTIFF}, ios::IO) = detectAspexTIFF(ios)

# ---------------------------------------------------------------------------------------------------#

const spectrumfiletypes = (ISOEMSA, BrukerSPX, ASPEXTIFF, BrukerPDZ)


function sniffspectrum(ios::IO)
    for st in spectrumfiletypes
        try
            seekstart(ios)
            if sniff(st, ios)
                seekstart(ios)
                return st
            end
        catch
            @info "Error sniffing filetype $st"
        end
    end
    @error "This object does not represent a known spectrum file type."
end

function sniffspectrum(filename::AbstractString)
    @assert isfile(filename) "$filename is not a file in sniffspectrum(...)."
    ext = lowercase(splitext(filename)[2])
    for st in spectrumfiletypes
        if ext in extensions(st)
            try
                res = open(filename, read = true) do ios
                    sniff(st, ios) ? st : nothing
                end
                if !isnothing(res)
                    return res
                end
            catch
                @info "Error sniffing filetype $st"
            end
        end
    end
    @error "This object does not represent a known spectrum file type."
end

loadspectrum(ios::IO) = loadspectrum(sniffspectrum(ios), ios)
loadspectrum(filename::AbstractString) = loadspectrum(sniffspectrum(filename), filename)
loadspectrum(filename::AbstractString, det::EDSDetector) =
    loadspectrum(sniffspectrum(filename), filename, det)
