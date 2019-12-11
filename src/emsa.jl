using FileIO

function isemsa(filename::AbstractString)
    open(filename,"r") do f
        for (lx, line) in enumerate(eachline(f))
        	res = split_emsa_header_item(line)
			if (res==nothing) || ((lx==1) && (res[1]!="FORMAT" || uppercase(res[2])!="EMSA/MAS SPECTRAL DATA FILE"))
                return false
			end
            if (res==nothing) || ((lx==2) && (res[1]!="VERSION" || res[2]!="1.0"))
                return false
            end
			if (lx>=2)
				return true
			end
		end
	end
end

function split_emsa_header_item(line::AbstractString)
    p = findfirst(":",line)
    if p ≠ nothing
        tmp = uppercase(strip(SubString(line,2:p.start-1)))
        pp = findfirst("-",tmp)
        if pp ≠ nothing
            key, mod = strip(SubString(tmp,1:pp.start-1)), strip(SubString(tmp,pp.stop+1))
        else
            key, mod = tmp, nothing
        end
        value = strip(SubString(line,p.stop+1))
        # println(key," = ",value)
        return (key, value, mod)
    else
        return nothing
    end
end

function parsecoating(value::AbstractString)::Union{Film,Missing}
    i=findfirst(" nm of ",value)
	if !isnothing(i)
		try
			thk = parse(Float64, value[1:i.start-1])*1.0e-7 # cm
			mat = parsedtsa2comp(value[i.stop+1:end])
			return Film(mat,thk)
		catch err
			@warn "Error parsing $(value) as a coating $(value)"
		end
	end
    return missing
end

"""
    readEMSA(filename::AbstractString, T::Type{<:Real}=Float64)::Spectrum

Read an ISO/EMSA format spectrum from a disk file at the specified path.
T is the type of the channel data elements.
"""
function readEMSA(filename::AbstractString, T::Type{<:Real}=Float64)::Spectrum
	open(filename) do f
		return readEMSA(f, T)
	end
end

"""
    readEMSA(ios::IOStream, T::Type{<:Real}=Float64)::Spectrum

Read an ISO/EMSA format spectrum from a disk file at the specified path.
T is the type of the channel data elements.
"""
function readEMSA(f::IOStream, T::Type{<:Real}=Float64)::Spectrum
    energy, counts = LinearEnergyScale(0.0,10.0), T[]
    props = Dict{Symbol,Any}()
    props[:Filename]=filename
    inData, xpcscale = 0, 1.0
    stgpos = Dict{Symbol,Float64}()
    for (lx, line) in enumerate(eachline(f))
        if (lx ≤ 2)
            res = split_emsa_header_item(line)
            if (res==nothing) || ((lx==1) && (res[1]!="FORMAT" || uppercase(res[2])!="EMSA/MAS SPECTRAL DATA FILE"))
                error("This file does not have the correct header to be an EMSA/MAS spectrum file.")
            end
            if (res==nothing) || ((lx==2) && (res[1]!="VERSION" || res[2]!="1.0"))
                error("This file is not a VERSION=1.0 EMSA/MAS spectrum file.")
            end
        elseif inData>0
            if startswith(line, "#ENDOFDATA")
                break
            else
                for cc in split(line, ",")
                    num = match(r"[-+]?[0-9]*\.?[0-9]+", cc)
                    if num ≠ nothing
                        append!(counts, parse(T, num.match))
                    end
                    inData=inData+1
                end
            end
        elseif startswith(line,"#")
            res = split_emsa_header_item(line)
            if res ≠ nothing
                key, value, mod = res
                if key == "SPECTRUM"
                    inData = 1
                elseif key == "BEAMKV"
                    props[:BeamEnergy]=1000.0*parse(Float64,value) # in eV
                elseif key == "XPERCHAN"
					xperch=parse(Float64,value)
					xpcscale = isnothing(mod) ? (xperch < 0.1 ? 1000.0 : 1.0) :  # Guess likely eV or keV
									(isequal(mod,"KEV") ? 1000.0 : 1.0)
                    energy = LinearEnergyScale(xpcscale*energy.offset, xpcscale*xperch)
                elseif key == "LIVETIME"
                    props[:LiveTime]=parse(Float64,value)
                elseif key == "REALTIME"
                    props[:RealTime]=parse(Float64,value)
                elseif key == "PROBECUR"
                    props[:ProbeCurrent]=parse(Float64,value)
                elseif key == "OFFSET"
                    energy = LinearEnergyScale(xpcscale*parse(Float64,value), energy.width)
                elseif key == "TITLE"
                    props[:Name]=value
                elseif key == "OWNER"
                    props[:Owner]=value
                elseif key == "ELEVANGLE"
                    props[:TakeOffAngle] = (props[:Elevation] = deg2rad(parse(Float64, value)))
                elseif key == "XPOSITION"
                    stgpos[:X] = parse(Float64, value)
                elseif key == "YPOSITION"
                    stgpos[:Y] = parse(Float64, value)
                elseif key == "ZPOSITION"
                    stgpos[:Z] = parse(Float64, value)
                elseif key == "COMMENT"
                    prev = getindex(props,:Comment,missing)
                    props[:Comment] = prev ≠ missing ? prev*"\n"*value : value
                elseif key == "#D2STDCMP"
                    props[:Composition] = parsedtsa2comp(value)
                elseif key == "#CONDCOATING"
                    props[:Coating] = parsecoating(value)
                elseif key =="#RTDET" # Bruker
					props[:DetectorModel] = value
				elseif key == "#TDEADLYR" # Bruker
					@assert isnothing(mod) || (!isequal(mod,"CM")) "Unexpected scale -$mod"
					props[:DeadLayer] = parse(Float64, value)
				elseif key == "#FANO" # Bruker
					props[:Fano] = parse(Float64, value)
				elseif key == "#MNFWHM" # Bruker
					sc = (isnothing(mod) ? 1000.0 : (isequal(mod,"KEV") ? 1000.0 : 1.0))
					props[:FWHMMnKa] = sc*parse(Float64, value)
				elseif key == "#IDENT" # Bruker
					elm = nothing
					try
						props[:XRFAnode] = parse(Element, value)
					catch
						# Just ignore it
					end
				end
            end
        end
    end
    if startswith( props[:Name], "Bruker")
        if haskey(props, :Composition)
            props[:Name]=name(props[:Composition])
        else
			path = splitpath(filename)
            props[:Name]=path[end]
        end
    end
	if endswith(props[:Name],".txt") || endswith(props[:Name],".msa")
		props[:Name] = props[:Name][1:end-4]
	end
    props[:StagePosition] = stgpos
    return Spectrum(energy, counts, props)
end

"""
    readEMSA(filename::AbstractString, det::Detector, force::Bool=false)::Spectrum

Read an EMSA file and apply the specified detector.  If force is false and the detector and
read calibration don't match then the function errors.
"""
function readEMSA(filename::AbstractString, det::Detector, force::Bool=false, T::Type{<:Real}=Float64)::Spectrum
	spec = readEMSA(filename, T)
	if !force
		if abs(channel(energy(1, det),spec)-1)>1
			error("Spectrum and detector calibrations don't match - Offset.")
		end
		test = det.channelcount÷2
		if abs(channel(energy(test,det),spec)-test)>1
			error("Spectrum and detector calibrations don't match - Gain.")
		end
	end
	return apply(spec, det)
end

const ISO_EMSA = format"ISO/EMSA"

function FileIO.load(file::File{ISO_EMSA})
    open(filename) do ios
        return load(Stream(ISO_EMSA,ios))
    end
end


function FileIO.load(ios::Stream{ISO_EMSA})
	return readEMSA(ios.io,Float64)
end

function FileIO.save(f::Stream{ISO_EMSA}, data)
    @error "Saving to ISO/EMSA streams is not implemented"
end

function FileIO.save(f::File{ISO_EMSA}, data)
    @error "Saving to ISO/EMSA files is not implemented"
end

FileIO.add_format(ISO_EMSA, isemsa, [".msa", ".emsa"])
