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

function writeEMSA(filename::AbstractString, spec::Spectrum)
	open(filename,"w") do ios
		writeEMSA(ios, spec)
	end
end

function writeEMSA(io::IOStream, spec::Spectrum)
	forceascii(ss) = map(c-> Int(c)&0x80==0 ? c : '?', ss)
	# "#FORMAT      : EMSA/MAS Spectral Data File"
	writeline(io, tag, data, extra="") =
		println(io,"#$(tag)$(repeat(' ',12-(length(tag)+length(extra))))$(extra): $(forceascii(data))")
	writeline(io, "FORMAT","EMSA/MAS Spectral Data File")
	writeline(io, "VERSION","1.0")
	writeline(io, "TITLE",spec[:Name])
	if haskey(spec,:AcquisitionTime)
		dt = spec[:AcquisitionTime]
		writeline(io, "DATE",uppercase(Dates.format(dt, "d-m-yyyy"))) # 25-FEB-2022
		writeline(io, "TIME",Times.format(dt, "HH:MM")) # 15:32
	end
	if haskey(spec,:Owner)
		writeline(io, "OWNER","")
	end
	writeline(io, "NPOINTS","$(length(spec))")
	writeline(io, "NCOLUMNS","1")
	writeline(io, "XUNITS","eV")
	writeline(io, "YUNITS","counts")
	if spec.energy isa LinearEnergyScale
		writeline(io, "DATATYPE","Y")
		writeline(io, "XPERCHAN","$(spec.energy.width)")
		writeline(io, "OFFSET", "$(spec.energy.offset)")
	else
		writeline(io, "DATATYPE","XY")
		dx = (energy(length(spec),spec)-energy(1,spec))/(length(spec)-1.0) # Mean width
		writeline(io, "XPERCHAN","$(dx)")
		writeline(io, "OFFSET", "$(energy(spec,1))") # energy of first channel
	end
	writeline(io, "SIGNALTYPE","EDS")
	if haskey(spec,:LiveTime)
		writeline(io, "LIVETIME","$(spec[:LiveTime])")
	end
	if haskey(spec,:RealTime)
		writeline(io, "REALTIME","$(spec[:RealTime])")
	end
	if haskey(spec,:BeamEnergy)
		writeline(io, "BEAMKV","$(0.001*spec[:BeamEnergy])")
	end
	if haskey(spec,:ProbeCurrent)
		writeline(io,"PROBECUR","$(spec[:ProbeCurrent])")
	end
	if haskey(spec,:TakeOffAngle)
		writeline(io, "ELEVANGLE","$(rad2deg(spec[:TakeOffAngle]))")
	end
	if haskey(spec,:Azimuthal)
		writeline(io, "AZIMANGLE","$(rad2deg(spec[:Azimuthal]))")
	end
	if haskey(spec,:StagePosition)
		sp=spec[:StagePosition]
		for (tag, sym) in ( ("XPOSITION", :X ), ("YPOSITION", :Y ), ("ZPOSITION", :Z ) )
			if haskey(sp, sym)
				writeline(io, tag, "$(10.0*sp[sym])","mm")
			end
		end
	end
	if haskey(spec,:Composition)
		writeline(io, "#D2STDCMP",todtsa2comp(spec[:Composition]))
	end
	if haskey(spec,:Specimen)
		writeline(io,"#SPECIMEN","$(spec[:Specimen])")
	end
	if haskey(spec,:Coating)
		cc = props[:Coating]
		writeline(io, "#CONDCOATING", "$(1.0e7*cc.thickness) nm of $(todtsa2comp(cc.material))")
	end
	if haskey(spec,:DetectorModel)
		writeline(io, "#RTDET",spec[:DetectorModel]) # Bruker
	end
	if haskey(spec,:DeadLayer)
		writeline(io, "#TDEADLYR", spec[:DeadLayer])  # Bruker
	end
	if haskey(spec,:Fano)
		writeline(io, "#FANO",spec[:Fano]) # Bruker
	end
	if haskey(spec, :FWHMMnKa)
		writeline(io, "#MNFWHM", spec[:FWHMMnKa], "EV") # Bruker
	end
	if haskey(spec, :XRFAnode)
		writeline(io, "#IDENT", spec[:XRFAnode]) # Bruker
	end
	if spec.energy isa LinearEnergyScale
		writeline(io, "SPECTRUM","Spectrum data follows in counts")
		for i in eachindex(spec)
			println(io,spec.counts[i])
		end
	else
		writeline(io, "SPECTRUM","Spectrum data follows in energy, counts")
		for i in eachindex(spec)
			println(io,"$(energy(i,spec)), $(spec.counts[i])")
		end
	end
	writeline(io, "#ENDOFDATA","")
end

const ISO_EMSA = format"ISO/EMSA"

function load(file::File{ISO_EMSA})
    open(filename) do ios
        return load(Stream(ISO_EMSA,ios))
    end
end

function load(ios::Stream{ISO_EMSA})
	return readEMSA(ios.io,Float64)
end

function save(f::Stream{ISO_EMSA}, spec::Spectrum)
    writeEMSA(f.io, spec)
end

function save(f::File{ISO_EMSA}, data)
	open(f) do ios
		writeEMSA(ios)
	end
end

FileIO.add_format(ISO_EMSA, isemsa, [".msa", ".emsa"])
