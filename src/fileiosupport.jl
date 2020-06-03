# This file defines a mechanism for auto-sniffing, loading and saving spectrum files

abstract type SpectrumFileType end

function savespectrum(ty::Type{<:SpectrumFileType}, fio, data)
    @error "Saving to is not implemented for $ty. Probably never will be."
end

loadspectrum(ty::Type{<:SpectrumFileType}, filename::AbstractString) =
    return open(filename, read=true) do ios
        loadspectrum(ty, ios)
    end

# ---------------------------------------------------------------------------------------------------#

struct BrukerSPX <: SpectrumFileType end

loadspectrum(::Type{BrukerSPX}, ios::IO) = readbrukerspx(ios)
extensions(::Type{BrukerSPX}) = ( ".spx", )


sniff(::Type{BrukerSPX}, ios::IO) = detectbrukerspx(ios)

# ---------------------------------------------------------------------------------------------------#

struct BrukerPDZ <: SpectrumFileType end

loadspectrum(::Type{BrukerPDZ}, ios::IO) = readbrukerpdz(ios)
extensions(::Type{BrukerPDZ}) = ( ".pdz", )

sniff(::Type{BrukerPDZ}, ios::IO) = detectbrukerpdz(ios)

# ---------------------------------------------------------------------------------------------------#

struct ISOEMSA <: SpectrumFileType end

loadspectrum(::Type{ISOEMSA}, ios::IO) = readEMSA(ios, Float64)
savespectrum(::Type{ISOEMSA}, ios::IO, spec::Spectrum) = writeEMSA(ios, spec)
savespectrum(::Type{ISOEMSA}, fn::AbstractString, spec::Spectrum) =
    open(fn,write=true) do ios
        writeEMSA(ios, spec)
    end
extensions(::Type{ISOEMSA}) = ( ".msa", ".emsa", ".ems", ".txt", )

sniff(::Type{ISOEMSA}, ios::IO) = isemsa(ios)

# ---------------------------------------------------------------------------------------------------#

struct ASPEXTIFF <: SpectrumFileType end

loadspectrum(::Type{ASPEXTIFF}, ios::IO) = readAspexTIFF(ios, withImgs=true)
extensions(::Type{ASPEXTIFF}) = ( ".tif", )

sniff(::Type{ASPEXTIFF}, ios::IO) = detectAspexTIFF(ios)

# ---------------------------------------------------------------------------------------------------#

const spectrumfiletypes = ( BrukerSPX, BrukerPDZ, ISOEMSA, ASPEXTIFF )


function sniffspectrum( ios::IO )
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

function sniffspectrum( filename::AbstractString )
    ext = lowercase(splitext(filename)[2])
    for st in spectrumfiletypes
        if ext in extensions(st)
            try
                res = open(filename, read=true) do ios
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
loadspectrum(filename::AbstractString) = loadspectrum(sniffspectrum(filename),filename)
