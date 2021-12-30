import FourierTools

# Keeps track of the number of spectra in this session.

let spectrumIndex = 0
    global spectrumCounter() = (spectrumIndex += 1)
end

"""
    Spectrum{T<:Real} <: AbstractVector{T}

    Spectrum(energy::EnergyScale, data::Vector{<:Real}, props::Dict{Symbol,Any})

Construct a structure to hold spectrum data (energy scale, counts and metadata).

See [`NeXLSpectrum.EnergyScale`](@ref) or [`NeXLSpectrum.LinearEnergyScale`](@ref)

Example:

    julia> spec = Spectrum(LinearEnergyScale(0.0,10.0),
                     collect(1:1024),
                     Dict{Symbol,Any}(:BeamEnergy=>10.0e3, :LiveTime=>30.0))
                                                                                  ** 1024.0
                                                                           *********
                                                                   *****************
                                                            ************************
                                                     *******************************
                                              **************************************
                                       *********************************************
                                ****************************************************
                         ***********************************************************
                  ******************************************************************
           *************************************************************************
    ******************************************************************************** 10.0 keV
    Spectrum[3062][10.0 keV, Unknown, 525000.0 counts]


Spectra are usually loaded from disk using:

    s1 = loadspectrum(joinpath(path,"ADM-6005a_1.msa")) # Auto-detects the file type (supports ISO/EMSA, Bruker, ASPEX, ??? files)
      *                                                                                                  128860.0
      *
      *
      *
      *
      *
      *     *
      *     *
      *     *
      * ** **
      * ** *****        **  **
    **************************************************************************************************** 20.0 keV
    Spectrum{Float64}[ADM-6005a_1, -484.20818 + 5.01716⋅ch eV, 4096 ch, 20.0 keV, Unknown, 6.81e6 counts]

`Spectrum` implements indexing using various different mechanisms.  If `spec` is a `Spectrum` then

    spec[123] # will return the number of counts in channel 123
	spec[123:222] # will return a Vector of counts from channel 123:222
    spec[134.] # will return the number of counts in the channel at energy 134.0 eV
    spec[134.0:270.0] # will return a Vector of counts for channels with energies between 134.0 eV and 270.0 eV
    spec[n"Ca L3-M5"] # Counts in the channel at the Ca L3-M5 characteristic X-ray energy
    spec[:Comment] # will return the meta-data item named :Comment

Basic spectrum math is supported using the operators *, /, ÷, +, and -.  If `s1`, `s2` and `s3` are `Spectrum` objects
then:

    2*s1 # A Spectrum containing twice as many counts in each channel
    2.0*s1 # A Spectrum containing twice as many counts in each channel
    s1/2.0 # A spectrum containing half as many counts in each channel
    s1÷2 # A spectrum containing half as many counts in each channel (Only works for Spectrum{<:Integer})
    s1 + s2 # A Spectrum containing channel-by-channel sum of counts in s1 and s1 and common properties.
    s1 - s2 # The channel-by-channel difference

Spectrum metadata is identified by a `Symbol`. These `Symbol`s are used within NeXLSpectrum. You can create others 
to associate other data items with a `Spectrum`.

    :BeamEnergy    # In eV
	:Elevation     # In radians
	:TakeOffAngle  # In radians (Detector position)
    :Azimuthal     # In radians (Detector position)
	:WorkingDistance # In cm (not mm!!!!)
    :LiveTime      # In seconds
    :RealTime      # In seconds
    :DeadFraction  # Fractional
    :ProbeCurrent  # In nano-amps
    :Name          # A `String`
    :Owner         # A `String`
    :Sample        # Description of the sample
    :StagePosition # A Dict{Symbol,Real} with entries :X, :Y, :Z, :R, :T, B: in cm and degrees
    :Comment       # A `String`
    :Composition   # A `Material` (known composition, not measured)
    :Elements      # A collection of elements in the material like `[ n"Fe", n"Ca", n"Si" ]`
    :Detector      # A Detector like a BasicEDS or another EDSDetector
    :Filename      # Source filename
    :Coating       # A `Film` or `Film[]` (eg. 10 nm of C|Au etc.)
	:AcquisitionTime # Date and time of acquisition (`DateTime` struct)
	:Signature     # Dict{Element,Real} with the "particle signature"
    :SolidAngle    # Detector solid angle is steradians (area/dist²)
    :Detector      # A detector model property like `EDSDetector`

Spectrum Image items:

    :ImageRotation # Rotation of the primary scan direction from the :X axis towards the :Y axis in radians

Less common items:

	:ImageMag	   # Magnification (assuming a 3.5\" image) of the first image
	:ImageZoom     # Additional zoom for second image in a two image TIFF
	:Operator      # Analyst in ASPEX TIFF files
    :Image1, :Image2 ... # Images associated with the spectrum
    :BrukerThroughtput # Nominal throughtput setting on a Bruker detector
    :DetectorSerialNumber # EDS detector serial number
    :DetectorModel # Vendor model name
    :DetectorThickness # Thickness of detector active area
    :DeadLayerThickness # Thickness of Si dead layer on the entrance surface of the detector
    :Window        # Window construction details
    :DetectorSolidAngle # Collection solid angle of the X-ray detector
    :ChamberPressure # Vacuum presure in the sample chamber
    :ChamberAtmosphere # Nominally the composition of the residual gas in the chamber

XRF related items:

    :BeamEnergy  # Accelarating voltage within X-ray tube (eV)
    :XRFTubeAnode    # Element from which the X-ray tube is constructed
    :ProbeCurrent  # Electron current in the X-ray tube
    :XRFTubeIncidentAngle # Incident angle of electron beam in tube
    :XRFTubeTakeOffAngle # Take-off angle from tube
    :XRFExcitationAngle # Angle of incidence of the X-ray beam on the sample
    :XRFDetectionAngle # Angle of the detector relative to the sample
    :XRFExcitationPathLength # Distance from X-ray source to sample
    :XRFDetectionPathLength # Distance from the sample to the X-ray detector
    :XRFSampleTilt    #  Additional tilt of the sample
    :XRFTubeWindow   # Construction of the tube window

Not all spectra will define all properties.  Algorithms can define the `NeXLCore.minproperties(ty::Type)` method to specify
which properties are required by an algorithm of `ty::Type`.  Then `hasminrequired` and `requiredbutmissing` methods
will determine whether a `Spectrum` or `Dict{Symbol,Any}` is suitable for an algorithm.
"""
struct Spectrum{T<:Real} <: AbstractVector{T}
    energy::EnergyScale
    counts::Vector{T}
    properties::Dict{Symbol,Any}
    hash::UInt  # Stays fixed as underlying data changes
    function Spectrum(energy::EnergyScale, data::AbstractVector{<:Real}, props::Dict{Symbol,<:Any})
        res = new{typeof(data[1])}(
            energy,
            data,
            convert(Dict{Symbol, Any}, props),
            hash(energy, hash(data, hash(props))),
        )
        res.properties[:Name] = get(res.properties, :Name, "Spectrum[$(spectrumCounter())]")
        res
    end
end

Base.hash(spec::Spectrum, h::UInt) = hash(spec.hash, h)
Base.isequal(spec1::Spectrum, spec2::Spectrum) =
    (spec1===spec2) || (
    (spec1.hash == spec2.hash) &&
    isequal(spec1.energy, spec2.energy) &&
    isequal(spec1.properties, spec2.properties) &&
    isequal(spec1.counts, spec2.counts))
Base.isless(s1::Spectrum, s2::Spectrum) =
    isequal(s1[:Name], s2[:Name]) ? isless(s1.hash, s2.hash) : isless(s1[:Name], s2[:Name])
# Make it act like an AbstractVector
Base.eltype(::Spectrum{T}) where { T <: Real } = T
Base.length(spec::Spectrum) = length(spec.counts)
Base.ndims(spec::Spectrum) = 1
Base.size(spec::Spectrum) = size(spec.counts)
Base.axes(spec::Spectrum) = Base.axes(spec.counts)
Base.eachindex(spec::Spectrum) = eachindex(spec.counts)
Base.strides(spec::Spectrum) = strides(spec.counts)
Base.firstindex(spec::Spectrum) = firstindex(spec.counts)
Base.lastindex(spec::Spectrum) = lastindex(spec.counts)
Base.IndexStyle(::Type{<:Spectrum{T}}) where {T <: Real} = IndexStyle(Vector{T}) 
Base.IteratorSize(::Type{<:Spectrum{T}}) where {T <: Real} =
    Base.IteratorSize(Vector{T})
Base.IteratorEltype(::Type{<:Spectrum{T}}) where {T <: Real} =
    Base.IteratorEltype(Vector{T})
Base.iterate(spec::Spectrum, args...) = iterate(spec.counts, args...)
# Integer indices
Base.@propagate_inbounds Base.getindex(spec::Spectrum, I...) =
    getindex(spec.counts, I...)
Base.get(spec::Spectrum, idx::Int, def = zero(eltype(spec.counts))) =
    get(spec.counts, idx, def)
# AbstractFloat indices
Base.getindex(spec::Spectrum, energy::AbstractFloat) = getindex(spec.counts, channel(energy, spec))
Base.getindex(spec::Spectrum, sr::StepRangeLen{<:AbstractFloat}) =
    spec.counts[channel(sr[1], spec):channel(sr[end], spec)]
# CharXRay indices
Base.getindex(spec::Spectrum, cxr::CharXRay) = spec.counts[channel(energy(cxr), spec)]
# Symbol indices (ie Properties)
Base.get(spec::Spectrum, sym::Symbol, def::Any = missing) = get(spec.properties, sym, def)
Base.getindex(spec::Spectrum, sym::Symbol)::Any = getindex(spec.properties, sym)
# setindex!(...) for channel data
Base.setindex!(spec::Spectrum{T}, val::Real, idx...) where {T <: Real} =
    setindex!(spec.counts, convert(T, val), idx...)
# setindex!(...) for spectrum properties
Base.setindex!(spec::Spectrum, val::Real, sym::Symbol) = 
    setindex!(spec.properties, val, sym)
Base.setindex!(spec::Spectrum, val::Any, sym::Symbol) = 
    setindex!(spec.properties, val, sym)
    
Base.copy(spec::Spectrum) = Spectrum(spec.energy, copy(spec.counts), copy(spec.properties))
Base.merge!(spec::Spectrum, props::Dict{Symbol,Any}) = merge!(spec.properties, props)
properties(spec::Spectrum) = spec.properties
NeXLCore.name(spec::Spectrum) = spec[:Name]
Base.convert(::Type{Spectrum{T}}, spec::Spectrum{T}) where { T <: Real } = spec
function Base.convert(::Type{Spectrum{T}}, spec::Spectrum{U})::Spectrum{T} where {T <: Real, U <: Real }
    return Spectrum(spec.energy, T.(spec.counts), copy(spec.properties))
end
function Base.similar(spec::Spectrum{T}, ::Type{U}) where { T <: Real, U <: Real }
    return Spectrum(spec.energy, similar(spec.counts, U), copy(spec.properties))
end

# Spectrum math simply performs the operation on the channel data but doesn't change any properties except the name.
Base.:*(a::Real, s::Spectrum) = property!(Spectrum(s.energy, a*s.counts, copy(s.properties)), :Name, "$a⋅$(s[:Name])")
Base.:*(s::Spectrum, a::Real) = property!(Spectrum(s.energy, a*s.counts, copy(s.properties)), :Name, "$a⋅$(s[:Name])")
Base.:/(s::Spectrum, a::Real) = property!(Spectrum(s.energy, s.counts/a, copy(s.properties)), :Name, "$(s[:Name])/$a")
Base.:÷(s::Spectrum{T}, a::Integer) where {T <: Integer} = property!(Spectrum(s.energy, s.counts .÷ a, copy(s.properties)), :Name, "$(s[:Name])÷$a")
# Operations with two Spectrum arguments combine the common properties and rename.
function Base.:+(s1::Spectrum, s2::Spectrum) 
    @assert isequal(s1.energy,s2.energy) "The energy calibration must be the same to add two spectra."
    r = Base.OneTo(min(length(s1.counts),length(s2.counts)))
    property!(Spectrum(s1.energy, view(s1.counts,r) + view(s2.counts,r), commonproperties(s1,s2)), :Name, "$(s1[:Name]) + $(s2[:Name])")
end
function Base.:-(s1::Spectrum, s2::Spectrum) 
    @assert isequal(s1.energy,s2.energy) "The energy calibration must be the same to subtract one spectrum from another."
    r = Base.OneTo(min(length(s1.counts),length(s2.counts)))
    property!(Spectrum(s1.energy, view(s1.counts,r) - view(s2.counts,r), commonproperties(s1,s2)), :Name, "$(s1[:Name]) - $(s2[:Name])")
end

"""
    offset(s::Spectrum, dcounts::Real)

Returns a `Spectrum` like `s` but with dcounts added to each channel.
"""
offset(s::Spectrum, dcounts::Real) = property!(Spectrum(s.energy, s.counts .+ dcounts, copy(s.properties)), :Name, "$(s[:Name]) offset by $dcounts")



"""
    property!(spec::Spectrum, sym::Symbol, val::Any)

Useful to broadcast properties over many spectra.
"""
function property!(spec::Spectrum, sym::Symbol, val::Any)
    spec[sym]=val
    spec
end

"""
    rangeofenergies(ch::Integer, spec::Spectrum)

Returns the low and high energy extremes for the channels `ch`.
"""
rangeofenergies(ch::Integer, spec::Spectrum) =
    (energy(ch, spec.energy), energy(ch + 1, spec.energy))

"""
    NeXLCore.hasminrequired(ty::Type, spec::Spectrum)

Does `spec` have the necessary property items for the algorithm `ty`?
"""
NeXLCore.hasminrequired(ty::Type, spec::Spectrum) = hasminrequired(ty, spec.properties)

"""
    NeXLCore.requiredbutmissing(ty::Type, spec::Spectrum)

Which properties are missing from `spec` but are required for the algorithm `ty`?
"""
NeXLCore.requiredbutmissing(ty::Type, spec::Spectrum) =
    requiredbutmissing(ty, spec.properties)

maxproperty(specs::AbstractArray{<:Spectrum}, prop::Symbol) =
    maximum(spec -> spec[prop], specs)
minproperty(specs::AbstractArray{<:Spectrum}, prop::Symbol) =
    minimum(spec -> spec[prop], specs)
sameproperty(specs::AbstractArray{<:Spectrum}, prop::Symbol) =
    all(spec -> spec[prop] == specs[1][prop], specs) ? #
    specs[1][prop] : #
    error("The property $prop is not equal for all these spectra.")

function Base.show(io::IO, ::MIME"text/plain", spec::Spectrum{<:Real})
    comp = haskey(spec, :Composition) ? name(spec[:Composition]) : "Unknown"
    e0 =
        haskey(spec, :BeamEnergy) ? #
        "$(round(spec[:BeamEnergy]/1000.0,sigdigits=3))" :  #
        "Unknown"
    cnts = "$(round(sum(spec.counts),sigdigits=3))"
    if !(get(io, :compact, false) || haskey(io, :SHOWN_SET))
        textplot(io, spec, size = (12, 100))
        print(io, "\n")
    end
    print(
        io,
        "Spectrum{$(eltype(spec.counts))}[$(spec[:Name]), $(spec.energy), $(length(spec.counts)) ch, $e0 keV, $comp, $cnts counts]",
    )
end

Base.show(io::IO, spec::Spectrum) = show(io, "text/plain", spec)

function textplot(io::IO, spec::Spectrum; size = (16, 80))
    (rows, cols) = size
    # how much to plot
    e0_eV = haskey(spec, :BeamEnergy) ? min(spec[:BeamEnergy], energy(length(spec), spec)) : energy(length(spec), spec)
    maxCh = min(channel(convert(Float64, e0_eV), spec), length(spec))
    step, max = maxCh ÷ cols, maximum(spec.counts)
    maxes = [rows * (maximum(spec.counts[(i-1)*step+1:i*step]) / max) for i in 1:cols]
    for r in 1:rows
        print(io, join(map(i -> r ≥ rows - maxes[i] ? '*' : ' ', 1:cols)))
        if r == 1
            println(io, " $(round(max,digits=0))")
        elseif r == rows
            print(io, " $(round(0.001*e0_eV,digits=2)) keV")
        else
            println(io)
        end
    end
end

textplot(spec::Spectrum; size = (16, 80)) = textplot(stdout, spec, size = size)

"""
    NeXLUncertainties.asa(::Type{DataFrame}, spec::Spectrum; properties::Bool = false)

Converts the spectrum energy and counts data into a DataFrame.
"""
NeXLUncertainties.asa(
    ::Type{DataFrame},
    spec::Spectrum;
    properties::Bool = false,
)::DataFrame =
    properties ?
    DataFrame(
        Keys = [keys(spec.properties)...],
        Values = repr.([values(spec.properties)...]),
    ) : #
    DataFrame(E = energyscale(spec), I = counts(spec))

"""
    apply(spec::Spectrum, det::EDSDetector)::Spectrum

Applies the specified detector to this spectrum by ensuring that the energy scales match and spec[:Detector] = det.
Creates a copy of the original spectrum unless the detector is already det.
"""

function apply(spec::Spectrum, det::EDSDetector)::Spectrum
    if (haskey(spec, :Detector)) && (spec[:Detector] == det)
        return spec
    else
        props = copy(spec.properties)
        props[:Detector] = det
        return Spectrum(det.scale, spec.counts, props)
    end
end

"""
	matching(spec::Spectrum, resMnKa::Float64, lld::Int=1, minByFam::Dict{Shell,Element} = Dict())::BasicEDS

Build an EDSDetector to match the channel count and energy scale in this spectrum.
"""
matching(spec::Spectrum, resMnKa::Float64, lld::Int = 1)::BasicEDS =
    BasicEDS(length(spec), spec.energy, MnKaResolution(resMnKa), lld)

matching(
    spec::Spectrum,
    resMnKa::Float64,
    lld::Int,
    minByFam::Dict{Shell,Element},
)::BasicEDS = BasicEDS(length(spec), spec.energy, MnKaResolution(resMnKa), lld, minByFam)
"""
    matches(spec::Spectrum, det::Detector, tol::Float64 = 1.0)::Bool
    matches(spec1::Spectrum, spec2::Spectrum, tol::Float64 = 1.0)::Bool

Does the calibration of the Spectrum (approximately) match the calibration of the Detector
or the other Spectrum?
"""
matches(spec::Spectrum, det::Detector, tol::Float64 = 1.0)::Bool = false
function matches(spec::Spectrum, det::EDSDetector, tol::Float64 = 1.0)::Bool
    return abs(energy(1, spec) - energy(1, det)) < tol * channelwidth(1, det) &&
           abs(energy(length(spec), spec) - energy(length(spec), det)) <
           tol * channelwidth(length(spec), det)
end
function matches(spec1::Spectrum, spec2::Spectrum, tol::Float64 = 1.0)::Bool
    return abs(energy(1, spec1) - energy(1, spec2)) < tol * channelwidth(1, spec1) &&
           abs(energy(length(spec1), spec1) - energy(length(spec2), spec2)) <
           tol * channelwidth(length(spec1), spec1)
end

"""
    Base.keys(spec::Spectrum)

Return the defined properties as a set of `Symbol`s.
"""
Base.keys(spec::Spectrum) = keys(spec.properties)

"""
    Base.haskey(spec::Spectrum, sym::Symbol)

Is the specified key defined?
"""
Base.haskey(spec::Spectrum, sym::Symbol) = haskey(spec.properties, sym)

"""
    elms(spec::Spectrum, withcoating = false)::Set{Element}

Returns a set of the elements associated with this spectrum. `withcoating` determines whether the coating
elements are also added.
"""
function NeXLCore.elms(spec::Spectrum, withcoating = false)::Set{Element}
    res = Set{Element}()
    if haskey(spec, :Elements)
        union!(res, spec[:Elements])
    end
    if haskey(spec, :Composition)
        union!(res, keys(spec[:Composition]))
    end
    if haskey(spec, :Signature)
        union!(res, keys(spec[:Signature]))
    end
    if withcoating && haskey(spec, :Coating)
        union!(res, keys(material(spec[:Coating])))
    end
    return res
end

"""
    dose(spec::Spectrum, def=missing)
    dose(props::Dict{Symbol, Any}, def=missing)

The probe dose in nA⋅s and is equals spec[:LiveTime]⋅spec[:ProbeCurrent].
"""
function dose(props::Dict{Symbol, Any}, def::Union{Float64,Missing})::Union{Float64,Missing}
    res = get(props, :LiveTime, missing) * get(props, :ProbeCurrent, missing)
    return isequal(res, missing) ? def : res
end
dose(spec::Spectrum, def::Union{Float64,Missing})::Union{Float64,Missing} = dose(spec.properties, def)
function dose(props::Dict{Symbol, Any})
    res = get(props, :LiveTime, missing) * get(props, :ProbeCurrent, missing)
    ismissing(res) && @error "One or more of the properties necessary to calculate the dose (:ProbeCurrent and :LiveTime) is not available."
    return res
end
dose(spec::Spectrum) = dose(spec.properties)

"""
    NeXLCore.energy(ch::Int, spec::Spectrum)::Float64

What energy is associated with `ch` in `spec`?
"""
NeXLCore.energy(ch::Int, spec::Spectrum)::Float64 = energy(ch, spec.energy)

"""
    coordinate(spec::Spectrum)
    coordinate(hss::Hyperspectrum, idx::NTuple{N, <Integer})
    coordinate(hss::HyperSpectrum, ci::CartesianIndex)

Returns the stage coordinate corresponding most closely to the :StagePosition (compensated
for the image offset in a `HyperSpectrum`) Defaults to :X=>0, :Y=>0

Nominally coordinates are:  :X, :Y, :Z, :R, :T, :B

Use with: `image2stage(...)` and `stage2image(...)` to calculate the stage coordinate of image pixels.
"""
function coordinate(spec::Spectrum)
    return get(spec, :StagePosition, Dict{Symbol,Float64}(:X=>0.0, :Y=>0.0))
end

"""
    channelwidth(ch::Int, spec::Spectrum)::Float64
    channelwidth(ch::Int, spec::HyperSpectrum)::Float64

Returns the width of the `ch` channel in eV.

    channelwidth(spec::Spectrum)::Float64

Returns the mean channel width.
"""
channelwidth(ch::Int, spec::Spectrum) = energy(ch + 1, spec) - energy(ch, spec)

channelwidth(spec::Spectrum) = (energy(length(spec), spec) - energy(1, spec))/length(spec) 

channel(eV::Float64, spec::Spectrum)::Int = channel(eV, spec.energy)
channel(::Type{Float64}, eV::Float64, spec::Spectrum)::Float64 = channel(Float64, eV, spec.energy)

"""
    counts(spec::Spectrum, ::Type{T}, applyLLD=false)::Vector{T} where {T<:Number}
    counts(spec::Spectrum, channels::AbstractUnitRange{<:Integer}, ::Type{T}, applyLLD=false)::Vector{T} where {T<:Number}

Creates a copy of the spectrum counts data as the specified Number type. If the spectrum has a :Detector
property then the detector's lld (low-level discriminator) and applyLLD=true then the lld is applied to the result
by setting all channels less-than-or-equal to det.lld to zero. This method throws an error if the counts data
cannot be converted without loss to the type `T`.
"""
function counts(
    spec::Spectrum,
    ::Type{T} = Float64,
    applyLLD = false,
)::Array{T} where {T<:Number}
    res = T.(spec.counts)
    if applyLLD && haskey(spec, :Detector)
        fill!(view(res, 1:lld(spec[:Detector])), zero(T))
    end
    return res
end
function counts(
    spec::Spectrum,
    channels::AbstractRange{<:Integer},
    ::Type{T} = Float64,
    applyLLD = false,
)::Array{T} where {T<:Real}
    if (first(channels) >= 1) && (last(channels) <= length(spec))
        res = T.(spec.counts[channels])
    else
        res = zeros(T, length(channels))
        r = max(first(channels), 1):min(last(channels), length(spec))
        res[first(r)-first(channels)+1:last(r)-first(channels)+1] = spec.counts[r]
    end
    if applyLLD && haskey(spec, :Detector) && (lld(spec) <= first(channels))
        fill!(view(res, 1:lld(spec)-first(channels)+1), zero(T))
    end
    return res
end
function counts(
    spec::Spectrum,
    channels::AbstractVector{<:Integer},
    ::Type{T} = Float64,
    applyLLD = false,
)::Array{T} where {T<:Real}
    res = T.(spec.counts[channels])
    if applyLLD && haskey(spec, :Detector) && (lld(spec) <= first(channels))
        fill!(view(res, 1:lld(spec)-first(channels)+1), zero(T))
    end
    return res
end

"""
    lld(spec::Spectrum)

Gets the low-level discriminator associated with this spectrum if there is one.
"""
lld(spec::Spectrum) =
    haskey(spec.properties, :Detector) ? lld(spec.properties[:Detector]) : 1

"""
	Base.findmax(spec::Spectrum, chs::AbstractRange{<:Integer})
    Base.findmax(spec::Spectrum)

Returns the (maximum intensity, channel index) over the specified range of channels
"""
function Base.findmax(spec::Spectrum, chs::AbstractRange{<:Integer})
    max = findmax(spec.counts[chs])
    return (max[1] + first(chs) - 1, max[2])
end
function Base.findmax(spec::Spectrum)
    last = min(
        haskey(spec, :BeamEnergy) ? channel(spec[:BeamEnergy], spec) : length(spec.counts),
        length(spec.counts),
    )
    return findmax(spec.counts[lld(spec):last])
end

"""
    integrate(spec::Spectrum, channels::AbstractUnitRange{<:Integer})
    integrate(spec::Spectrum, energyRange::StepRangeLen{Float64})

Sums all the counts in the specified channels.  No background correction.

    integrate(spec::Spectrum, back1::AbstractUnitRange{<:Integer}, peak::AbstractUnitRange{<:Integer}, back2::AbstractUnitRange{<:Integer})::UncertainValue
    integrate(spec::Spectrum, back1::StepRangeLen{Float64}, peak::StepRangeLen{Float64}, back2::StepRangeLen{Float64})::UncertainValue

Sums all the counts in the specified channels with background correction using the background intervals.

    integrate(spec::Spectrum)

Total integral of all counts from the LLD to the beam energy

"""
integrate(spec::Spectrum, channels::AbstractUnitRange{<:Integer}) =
    sum(spec.counts[channels])

integrate(spec::Spectrum, energyRange::StepRangeLen{Float64}) =
    sum(spec.counts[channel(energyRange[1], spec):channel(energyRange[end], spec)])

function integrate(
    spec::Spectrum,
    back1::AbstractUnitRange{<:Integer},
    peak::AbstractUnitRange{<:Integer},
    back2::AbstractUnitRange{<:Integer},
)::UncertainValue
    bL, bH = mean(spec.counts[back1]), mean(spec.counts[back2])
    cL, cH, cP = mean(back1), mean(back2), mean(peak)
    iP, w, t = sum(spec.counts[peak]), length(peak), ((cP - cL) / (cH - cL))
    return uv(
        iP - w * (bL * (1.0 - t) + bH * t),
        sqrt(iP + (sqrt(bL) * w * (1.0 - t))^2 + (sqrt(bH) * w * t)^2),
    )
end

function integrate(
    spec::Spectrum,
    back1::StepRangeLen{Float64},
    peak::StepRangeLen{Float64},
    back2::StepRangeLen{Float64},
)::UncertainValue
    b1i = channel(back1[1], spec):channel(back1[end], spec)
    b2i = channel(back2[1], spec):channel(back2[end], spec)
    pi = channel(peak[1], spec):channel(peak[end], spec)
    integrate(spec, b1i, pi, b2i)
end

function integrate(spec::Spectrum)
    last = min(
        haskey(spec, :BeamEnergy) ? channel(spec[:BeamEnergy], spec) : length(spec.counts),
        length(spec.counts),
    )
    return integrate(spec, lld(spec):last)
end

"""
    kratio(unk::Spectrum, std::Spectrum, back1::AbstractUnitRange{<:Integer}, peak::AbstractUnitRange{<:Integer}, back2::AbstractUnitRange{<:Integer})::UncertainValue
    kratio(unk::Spectrum, std::Spectrum, back1::StepRangeLen{Float64}, peak::StepRangeLen{Float64}, back2::StepRangeLen{Float64})::UncertainValue

The k-ratio of unk relative to std corrected for dose.  Requires that `unk` and `std` have the properties
`:LiveTime` and `:ProbeCurrent` defined.
"""
function kratio(
    unk::Spectrum,
    std::Spectrum,
    back1::AbstractUnitRange{<:Integer},
    peak::AbstractUnitRange{<:Integer},
    back2::AbstractUnitRange{<:Integer},
)
    iunk, istd = integrate(unk, back1, peak, back2) / dose(unk),
    integrate(std, back1, peak, back2) / dose(std)
    f = NeXLUncertainties.value(iunk) / NeXLUncertainties.value(istd)
    return uv(
        f,
        abs(f) * sqrt(
            (σ(iunk) / NeXLUncertainties.value(iunk))^2 +
            (σ(istd) / NeXLUncertainties.value(istd))^2,
        ),
    )
end

kratio(
    unk::Spectrum,
    std::Spectrum,
    back1::StepRangeLen,
    peak::StepRangeLen,
    back2::StepRangeLen,
) = kratio(
    unk,
    std, #
    channel(back1[1], unk):channel(back1[end], unk), #
    channel(peak[1], unk):channel(peak[end], unk), #
    channel(back2[1], unk):channel(back2[end], unk),
)


"""
    energyscale(spec::Spectrum)
    energyscale(spec::Spectrum, chs::AbstractRange{<:Integer})

Returns an array with the bin-by-bin energies
"""
energyscale(spec::Spectrum) = energyscale(spec.energy, eachindex(spec))
energyscale(spec::Spectrum, chs::AbstractRange{<:Integer}) = energyscale(spec.energy, chs)

"""
    simpleEDS(spec::Spectrum, fwhmatmnka::Float64)

Build a `EDSDetector` object for this spectrum with the specified FWHM at Mn Kα.
"""
simpleEDS(spec::Spectrum, fwhmatmnka::Float64) =
    BasicEDS(length(spec), spec.energy, MnKaResolution(fwhmatmnka))


"""
    subsample(spec::Spectrum, frac::Float64)

Subsample the counts data in a spectrum according to a statistically valid algorithm.  Returns
`spec` if frac>=1.0.
"""
function subsample(spec::Spectrum, frac::Float64)::Spectrum
    @assert frac > 0.0 "frac must be larger than zero."
    @assert frac <= 1.0 "frac must be less than or equal to 1.0."
    @assert haskey(spec, :LiveTime) "Please specify a :LiveTime in subsample"
    ss(n, f) = n > 0 ? sum(rand() <= f ? 1 : 0 for _ in 1:n) : 0
    frac = max(0.0, min(frac, 1.0))
    if frac < 1.0
        props = deepcopy(spec.properties)
        props[:LiveTime] = frac * spec[:LiveTime] # Must have
        if haskey(spec, :RealTime)
            props[:RealTime] = frac * spec[:RealTime] # Might have
        end
        return Spectrum(spec.energy, map(n -> ss(floor(Int, n), frac), spec.counts), props)
    else
        @warn "Not actually subsampling the spectrum because $frac > 1.0."
        return spec
    end
end
"""
    subdivide(spec::Spectrum, n::Int)::Vector{Spectrum}

Splits the event data in one spectrum into n spectra by assigning each event
to a pseudo-random choice of one of the n result spectra.  Produces n spectra
that act as though the original spectrum was collected in n time intervals
of LiveTime/n.  This is quite slow because it needs to call rand() for each
count in the spectrum (not just each channel).
"""
function subdivide(spec::Spectrum{T}, n::Int)::Vector{Spectrum{T}} where {T <: Real}
    res = zeros(Int, n, length(spec.counts))
    for ch in eachindex(spec.counts)
        # Assign each event to one and only one detector
        si = rand(1:n, floor(Int, spec[ch]))
        for i in 1:n
            res[i, ch] = count(e -> e == i, si)
        end
    end
    return map(1:n) do i
        props = deepcopy(spec.properties)
        props[:Name] = "Sub[$(spec[:Name]),$(i) of $(n)]"
        props[:LiveTime] = spec[:LiveTime] / n # Must have
        if haskey(spec, :RealTime)
            props[:RealTime] = spec[:RealTime] / n # Might have
        end
        Spectrum(spec.energy, res[i, :], props)
    end
end


"""
    estimatebackground(data::AbstractArray{Float64}, channel::Int, width::Int=5, order::Int=2)

Returns the tangent to the a quadratic fit to the counts data centered at channel with width
"""
function estimatebackground(
    data::AbstractArray{Float64},
    channel::Int,
    width::Int = 5,
    order::Int = 2,
)::ImmutablePolynomial
    minCh, maxCh = max(1, channel - width), min(length(data), channel + width)
    if maxCh - minCh >= order
        fr = fit(
            ImmutablePolynomial,
            (minCh-channel):(maxCh-channel),
            data[minCh:maxCh],
            order,
        )
        return ImmutablePolynomial([fr(0), derivative(fr)(0)]) # Linear
    else
        return ImmutablePolynomial([mean(data[minCh:maxCh]), 0.0])
    end
end

"""
    modelBackground(spec::Spectrum, chs::AbstractUnitRange{<:Integer}, ash::AtomicSubShell)

spec: A spectrum containing a peak centered on chs
chs:  A range of channels containing a peak
ash:  The edge (as an AtomicSubShell)

A simple model for modeling the background under a characteristic x-ray peak. The model
fits a line to low and high energy background regions around firsr(chs) and last(chs). If
the low energy line extended out to the edge energy is larger than the high energy line
at the same energy, then a negative going edge is fit between the two. Otherwise a line
is fit between the low energy side and the high energy side. This model only works when
there are no peak interference over the range chs.


    modelBackground(spec::Spectrum, chs::AbstractUnitRange{<:Integer})
    modelBackground(spec::Spectrum, chs::AbstractUnitRange{<:Integer}, ash::AtomicSubShell)

`spec`: A spectrum containing a peak centered on chs
`chs`:  A range of channels containing a peak
`ash`:  The largest edge within the range of channels `chs` associated with the characteristic peak

A simple model for modeling the background under a characteristic x-ray peak. The model
fits a line between the  low and high energy background regions around first(chs) and last(chs).
This model only works when there are no peak interference over the range chs.
"""
function modelBackground(
    spec::Spectrum,
    chs::AbstractUnitRange{<:Integer},
    ash::AtomicSubShell,
)
    cnts, ec = counts(spec), channel(energy(ash), spec)
    bl, bh = estimatebackground(cnts, first(chs), 5), estimatebackground(cnts, last(chs), 5)
    # bh = ch-> mean(cnts[last(chs):min(length(cnts),last(chs)+5)])
    if (ec < last(chs)) &&
       (bl(ec - first(chs)) > bh(ec - last(chs))) &&
       (energy(ash) < 2.0e3)  && ( shell(ash.subshell)==KShell )
        res = zeros(Float64, length(chs))
        res[1:ec-first(chs)] .= (bl(y - 1) for y in 1:ec-first(chs))
        res[ec-first(chs)+1] = 0.5 * (bl(ec - first(chs) + 1) + bh(ec - last(chs)))
        res[ec-first(chs)+2:last(chs)-first(chs)+1] .=
            (bh(y - last(chs)) for y in ec+1:last(chs))
    else
        s = (bh(0) - bl(0)) / length(chs)
        back = ImmutablePolynomial([bl(0), s])
        res = back.(collect(0:length(chs)-1))
    end
    return res
end

function modelBackground(spec::Spectrum, chs::AbstractUnitRange{<:Integer})
    cxs = counts(spec)
    bl = estimatebackground(cxs, first(chs), 5)
    bh = estimatebackground(cxs, last(chs), 5)
    s = (bh(0) - bl(0)) / length(chs)
    back = ImmutablePolynomial([bl(0), s])
    return back.(collect(0:length(chs)-1))
end

"""
    extractcharacteristic(spec::Spectrum, lowBack::AbstractUnitRange{<:Integer}, highBack::AbstractUnitRange{<:Integer})::Vector{Float64}

Extract the characteristic intensity for the peak located within chs with an edge at ash.
"""
function extractcharacteristic(
    spec::Spectrum,
    chs::AbstractUnitRange{<:Integer},
    ash::AtomicSubShell,
)::Vector{Float64}
    return counts(spec, chs, Float64) - modelBackground(spec, chs, ash)
end


"""
    peak(spec::Spectrum, chs::AbstractUnitRange{<:Integer}, ash::AtomicSubShell)::Float64

Estimates the peak intensity for the characteristic X-ray in the specified range of channels.
"""
peak(spec::Spectrum, chs::AbstractUnitRange{<:Integer})::Float64 =
    return sum(counts(spec, chs, Float64)) - background(spec, chs)

"""
    background(spec::Spectrum, chs::AbstractUnitRange{<:Integer}, ash::AtomicSubShell)::Float64

Estimates the background intensity for the characteristic X-ray in the specified range of channels.
"""
background(spec::Spectrum, chs::AbstractUnitRange{<:Integer})::Float64 =
    sum(modelBackground(spec, chs))


"""
    peaktobackground(spec::Spectrum, chs::AbstractUnitRange{<:Integer}, ash::AtomicSubShell)::Float64

Estimates the peak-to-background ratio for the characteristic X-ray intensity in the specified range of channels
which encompass the specified AtomicSubShell.
"""
function peaktobackground(
    spec::Spectrum,
    chs::AbstractUnitRange{<:Integer},
    ash::AtomicSubShell,
)::Float64
    back = sum(modelBackground(spec, chs, ash))
    return (sum(counts(spec, chs, Float64)) - back) / back
end

"""
    estkratio(unk::Spectrum, std::Spectrum, chs::AbstractUnitRange{<:Integer})

Estimates the k-ratio from niave models of peak and background intensity.  Only works if the peak is not interfered.
"""
estkratio(unk::Spectrum, std::Spectrum, chs::AbstractUnitRange{<:Integer}) =
    peak(unk, chs) * dose(std) / (peak(std, chs) * dose(unk))


"""
    NeXLUncertainties.asa(::Type{DataFrame}, spec::AbstractArray{Spectrum})::DataFrame

Returns a DataFrame that summarizes the list of spectra.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    specs::AbstractArray{<:Spectrum},
)::DataFrame
    _asname(comp) = ismissing(comp) ? missing : name(comp)
    unf, unl, uns = Union{Float64,Missing}, Union{Film,Nothing}, Union{String,Missing}
    nme, e0, pc, lt, rt, coat, integ, comp =
        String[], unf[], unf[], unf[], unf[], unl[], Float64[], uns[]
    for spec in specs
        push!(nme, spec[:Name])
        push!(e0, get(spec, :BeamEnergy, missing))
        push!(pc, get(spec, :ProbeCurrent, missing))
        push!(lt, get(spec, :LiveTime, missing))
        push!(rt, get(spec, :RealTime, missing))
        push!(coat, get(spec, :Coating, nothing))
        push!(integ, integrate(spec))
        push!(comp, _asname(get(spec, :Composition, missing)))
    end
    return DataFrame(
        Name = nme,
        BeamEnergy = e0,
        ProbeCurrent = pc,
        LiveTime = lt,
        RealTime = rt,
        Coating = coat,
        Integral = integ,
        Material = comp,
    )
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, spec::AbstractDict{Element, Spectrum})::DataFrame

Returns a DataFrame that summarizes a dictionary of standard spectra.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    stds::AbstractDict{Element,Spectrum},
)::DataFrame
    _asname(comp) = ismissing(comp) ? missing : name(comp)
    unf, unl, uns = Union{Float64,Missing}, Union{Film,Nothing}, Union{String,Missing}
    elm, zs, mfs, nme, e0, pc, lt, rt, coat, integ, comp = String[],
    Int[],
    Float64[],
    String[],
    unf[],
    unf[],
    unf[],
    unf[],
    unl[],
    Float64[],
    uns[]
    for el in sort(collect(keys(stds)))
        spec = stds[el]
        push!(elm, el.symbol)
        push!(zs, el.number)
        push!(nme, spec[:Name])
        push!(mfs, ismissing(spec[:Composition]) ? missing : spec[:Composition][el])
        push!(e0, get(spec, :BeamEnergy, missing))
        push!(pc, get(spec, :ProbeCurrent, missing))
        push!(lt, get(spec, :LiveTime, missing))
        push!(rt, get(spec, :RealTime, missing))
        push!(coat, get(spec, :Coating, nothing))
        push!(integ, integrate(spec))
        push!(comp, _asname(get(spec, :Composition, missing)))
    end
    return DataFrame(
        Element = elm,
        Z = zs,
        Name = nme,
        Material = comp,
        MassFrac = mfs,
        BeamEnergy = e0,
        ProbeCurrent = pc,
        LiveTime = lt,
        RealTime = rt,
        Coating = coat,
        Integral = integ,
    )
end

"""
    details(io, spec::Spectrum)
    details(spec::Spectrum)::String

Outputs a description of the data in the spectrum.
"""
function details(spec::Spectrum)
    df = DataFrame(Property=String[], Value=String[])
    push!(df, ("Name", "$(spec[:Name])"))
    push!(df, ("Calibration", "$(spec.energy)"))
    push!(df, ("Beam energy", "$(get(spec, :BeamEnergy, missing)/1000.0) keV"))
    push!(df, ("Probe current", "$(get(spec, :ProbeCurrent, missing)) nA"))
    push!(df, ("Live time", "$(get(spec, :LiveTime, missing)) s"))
    push!(df, ("Coating", "$(get(spec,:Coating, nothing))"))
    push!(df, ("Detector", "$(get(spec, :Detector, missing))"))
    push!(df, ("Comment", "$(get(spec, :Comment, missing))"))
    push!(df, ("Integral", "$(integrate(spec)) counts"))
    comp = get(spec, :Composition, missing)
    if !ismissing(comp)
        push!(df, ("Composition", "$(comp)"))
        det = get(spec, :Detector, missing)
        if !ismissing(det)
            coating = get(spec, :Coating, missing)
            comp2 = Set(keys(comp))
            if !ismissing(coating)
                union!(comp2, keys(coating.material))
            end
            for elm1 in comp2
                for ext1 in extents(elm1, det, 1.0e-4)
                    intersects = []
                    for elm2 in comp2
                        if elm2 ≠ elm1
                            for ext2 in extents(elm2, det, 1.0e-4)
                                if length(intersect(ext1, ext2)) > 0
                                    push!(intersects, "$(elm2.symbol)[$(ext2)]")
                                end
                            end
                        end
                    end
                    if length(intersects) > 0
                        push!(df, ("ROI $(elm1.symbol)[$(ext1)]","Intersects $(join(intersects,", "))"))
                    else
                        p, b = peak(spec, ext1), background(spec, ext1)
                        σ = p / sqrt(b)
                        push!(df, ("ROI $(elm1.symbol)[$(ext1)]", "$(round(Int,p)) counts over $(round(Int,b)) counts - σ = $(round(Int,σ))"))
                    end
                end
            end
        end
    end
    return df
end

details(io::IO, spec::Spectrum) = print(io, details(spec))

"""
    commonproperties(specs::AbstractArray{Spectrum})
    commonproperties(props1::Dict{Symbol,Any}, props2::Dict{Symbol,Any})

Return the properties that are held in common by all the spectra.
"""
function commonproperties(props1::Dict{Symbol,Any}, props2::Dict{Symbol,Any})
    res = Dict{Symbol,Any}()
    for (key, sp1) in props1
        if isequal(sp1, get(props2, key, missing))
            res[key] = sp1
        end
    end
    return res
end
function commonproperties(specs::Vararg{<:Spectrum})
    props = specs[1].properties
    for sp2 in specs[2:end]
        props = commonproperties(props, sp2.properties)
    end
    return props
end
commonproperties(specs::AbstractArray{<:Spectrum}) = commonproperties(specs...)

"""
    maxspectrum(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where T<:Real

Compute the max-pixel spectrum for the specified spectra.
"""
function maxspectrum(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where {T<:Real}
    res = zeros(T, maximum(length.(specs)))
    for spec in specs
        for i in eachindex(spec)
            res[i] = max(res[i], spec[i])
        end
    end
    props = commonproperties(specs)
    props[:Name] = "MaxSpectrum[$(length(specs)) spectra]"
    return Spectrum(spec1.energy, res, props)
end

"""
    minspectrum(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where T<:Real

Compute the min-pixel spectrum for the specified spectra.
"""
function minspectrum(specs::AbstractArray{Spectrum{T}})::Spectrum{T} where {T<:Real}
    res = zeros(T, maximum(length.(specs)))
    fill!(res, typemax(T))
    for spec in specs
        for i in eachindex(spec)
            res[i] = min(res[i], spec[i])
        end
    end
    props = commonproperties(specs)
    props[:Name] = "MinSpectrum[$(length(specs)) spectra]"
    return Spectrum(spec1.energy, res, props)
end


"""
    Base.sum(specs::AbstractArray{Spectrum{T}}; restype::Type{<:Real}=T, applylld=false, name=Nothing|String)::Spectrum{T}

Computes the sum spectrum over an `AbstractArray{Spectrum}` where the :ProbeCurrent and :LiveTime will be maintained
in a way that maintains the sum of the individual doses.  This function assumes (but does not check) that the energy
scales are equivalent for all the spectra.  The resultant energy scale is the scale of the first spectrum.  Other than
:ProbeCurrent, :LiveTime and :RealTime, only the properties that the spectra hold in common will be maintained.

This function behaves differently from `reduce(+, specs)` which checks whether the energy scales match and fails if 
they don't.
"""
function Base.sum(
    specs::AbstractArray{Spectrum{T}},
    restype::Type{<:Real} = T;
    applylld = false,
    name = nothing,
)::Spectrum where {T<:Real}
    cxs = zeros(restype, maximum(length.(specs)))
    for spec in specs
        cxs .+= counts(spec, eachindex(cxs), restype, applylld)
    end
    props = commonproperties(specs)
    props[:ProbeCurrent] = mean(sp[:ProbeCurrent] for sp in specs)
    props[:LiveTime] = sum(dose.(specs)) / props[:ProbeCurrent]
    rt = sum(get(sp, :RealTime, NaN64) for sp in specs)
    if !isnan(rt)
        props[:RealTime] = rt
    end
    props[:Name] = something(name, "Sum[$(length(specs)) spectra]")
    return Spectrum(specs[1].energy, cxs, props)
end

"""
    χ²(s1::Spectrum{T}, s2::Spectrum{T}, chs)::T where {T<:Real}
    χ²(specs::AbstractArray{Spectrum{T}}, chs)::Matrix{T}

Computes the dose corrected χ² metric of similarity between `s1` and `s2`.
Or computes the matrix of χ² for the array of spectra.
"""
function χ²(s1::Spectrum{T}, s2::Spectrum{T}, chs)::T where {T<:Real}
    k1, k2 = 1.0 / dose(s1), 1.0 / dose(s2)
    d = sum((k1 * s1[ch] - k2 * s2[ch])^2 for ch in chs)
    s = sum(sqrt(k1^2 * max(one(T), s1[ch]) + k2^2 * max(one(T), s2[ch])) for ch in chs)
    return d / s
end

function χ²(specs::AbstractVector{Spectrum{T}}, chs)::Matrix{T} where {T<:Real}
    χ2s = zeros(Float64, (length(specs), length(specs)))
    for i in eachindex(specs), j in i+1:length(specs)
        χ2s[i, j] = χ²(specs[i], specs[j], chs)
        χ2s[j, i] = χ2s[i, j]
    end
    return χ2s
end

"""
    measure_dissimilarity(specs::AbstractArray{Spectrum{T}}, chs)::Vector{Float64}
    measure_dissimilarity(specs::AbstractArray{Spectrum}, minE::Float64=100.0)::Vector{Float64}
    measure_dissimilarity(specs::AbstractArray{Spectrum{T}}, det::Detector, mat::Material)::Vector{Float64}

Returns a vector of χ²s which measure how different the i-th `Spectrum` is from the other spectra.
The first version covers all the channels between minE and the nominal beam energy. The second version
only considers those channels representing peaks in a spectrum from the `Material` on the `Detector`.
The statistic returned is the χ² for each spectrum relative to the mean spectrum.  It should be nominally
a vector of 1's for spectra that differ only by count statistics.
"""
function measure_dissimilarity(
    specs::AbstractArray{<:Spectrum},
    chs,
)::Vector{Float64}
    s = sum(specs)
    return [χ²(spec, s, chs) for spec in specs]
end

function measure_dissimilarity(
    specs::AbstractArray{<:Spectrum},
    minE::Float64 = 100.0,
)::Vector{Float64}
    e0 = maximum(spec[:BeamEnergy] for spec in specs)
    chs =
        minimum(
            channel(minE, spec) for spec in specs
        ):maximum(channel(e0, spec) for spec in specs)
    return measure_dissimilarity(specs, chs)
end

function measure_dissimilarity(
    specs::AbstractArray{<:Spectrum},
    det::Detector,
    mat::Material,
)::Vector{Float64}
    function mrg(inp::Vector{UnitRange{Int}})::Vector{UnitRange{Int}}
        simp = sort(inp)
        st, res = simp[1], UnitRange{Int}[]
        for r in simp[2:end]
            if isempty(intersect(st, r))
                push!(res, st)
                st = r
            else
                st = first(st):max(last(r), last(st))
            end
        end
        push!(res, st)
        return res
    end
    e0 = maximum(spec[:BeamEnergy] for spec in specs)
    # Figure out the contiguous ROIs and the channels in the ROIs 
    rois = mrg(
        mapreduce(
            elm -> extents(characteristic(elm, alltransitions, 0.01, e0), det, 0.001),
            append!,
            keys(mat),
        ),
    )
    chs = mapreduce(collect, append!, rois)
    return measure_dissimilarity(specs, chs)
end


"""
    findsimilar(specs::AbstractArray{Spectrum{T}}; tol = 1.8, minspecs=3)::Vector{Spectrum{T}}

Filters a collection of spectra for the ones most similar to the average by
removing the least similar spectrum sequentially until all the remaining spectra are within
tol times the mean(χ²).  This is useful for finding which of a set of replicate spectra are
sufficiently similar to each other.
"""
function findsimilar(
    specs::AbstractArray{Spectrum{T}};
    tol = 1.8,
    minspecs = 3,
)::Vector{Spectrum{T}} where {T<:Real}
    if length(specs) > minspecs
        σs = measure_dissimilarity(specs)
        (fmσ, fmi) = findmax(σs)
        rem = filter(σ -> σ ≠ fmσ, σs)
        # Now perform this recursively until all are within tol or we hit minspecs
        if fmσ > tol * Statistics.mean(rem)
            return findsimilar(
                filter(s -> s != specs[fmi], specs),
                tol = tol,
                minspecs = minspecs,
            )
        end
    end
    return specs
end


# Duane-Hunt related functions

_duane_hunt_func(es, p) = map(ee-> p[1]*bremsstrahlung(Small1987, ee, p[2], n"H"), es)

function _duane_hunt_impl(spec::Spectrum)
    if channel(spec[:BeamEnergy], spec) < length(spec)
        cdata = counts(spec)
        este0, estCh = let # Sometimes the D-H can be much below `spec[:BeamEnergy]`
            # this handles pulse pile-up above the D-H
            lim = min(10.0, max(1.0, 1.0e-5 * maximum(cdata)))
            len = min(10, max(5, length(spec) ÷ 400))
            chlast = findlast(i -> mean(@view cdata[i-len:i]) > lim, len+1:length(cdata))
            energy(chlast, spec), chlast
        end
        xdata, ydata = let
            # Range of channels above/below `este0`
            chs = max(1,(9*estCh)÷10):min((21*estCh)÷20,length(spec))
            energyscale(spec, chs), cdata[chs]
        end
        esti0 = sum(ydata) / sum(_duane_hunt_func(xdata, (1.0, este0)))
        e0 = spec[:BeamEnergy]
        # Apply constraint functions to keep DH between 0.1 E0 and 1.1 E0.
        transform(v) = asin(2*(v-0.1*e0)/(1.1*e0-0.1*e0)-1)
        inv_transform(tv) = 0.1*e0+(sin(tv)+1)*(1.1*e0-0.1*e0)/2
        # Fit Small's continuum model to the data (ignore detector efficiency)
        model(xs, p) = _duane_hunt_func(xs, (p[1], inv_transform(p[2])))
        res = curve_fit(model, xdata, ydata, [esti0, transform(este0)])
        return (res.param[1], inv_transform(res.param[2]))
    else
        error("Unable to estimate the Duane-Hunt on $(spec[:Name]) because the counts data does not extend to $(spec[:BeamEnergy]) keV.")
    end
end

"""
    duane_hunt(spec::Spectrum)

Estimates the Duane-Hunt limit (the energy at which the continuum goes to effectively zero.)
"""
duane_hunt(spec::Spectrum)  = _duane_hunt_impl(spec)[2]

"""
    uv(spec::Spectrum, chs::AbstractRange{<:Integer}=eachindex(spec))::Vector{UncertainValue}

Converts the count's data in a spectrum into an Vector{UncertainValue} assuming
count statistics can be approximated by C ± √C. 
"""
function NeXLUncertainties.uv(spec::Spectrum, chs::AbstractRange{<:Integer}=eachindex(spec))::Vector{UncertainValue}
    c = counts(spec, chs)
    return uv.(c, sqrt.(max.(c,one(eltype(spec)))))
end

"""
    sigma(spec::Spectrum, specs::AbstractArray{<:Spectrum}, chs::AbstractRange{<:Integer})::Vector{Float64}

Computes on a channel-by-channel basis how much `spec` spectrum deviates from the mean of the
other spectra in `specs`.  The result is expressed in terms of the standard deviation expected
from count statistics alone.   Assuming `spec` varies only by count statistics we expect
the result values have a mean 0.0 and a standard deviation of 1.0. 
"""
function sigma(spec::Spectrum, specs::AbstractArray{<:Spectrum}, chs::AbstractRange{<:Integer})::Vector{Float64}
    function doseaverage(specs, chs) 
        t = [ uv(spec, chs)/dose(spec) for spec in specs ]
        return map(j->mean(collect(t[i][j] for i in eachindex(t))), eachindex(t[1]))
    end
    function delta(spec, specs, chs)
        minus(uv1,uv2) = uv(value(uv1)-value(uv2),sqrt(variance(uv1)+variance(uv2)))
        return minus.(uv(spec, chs), dose(spec)*doseaverage(filter(s->s!=spec, specs), chs))
    end
    return map(v->value(v)/σ(v), delta(spec, specs, chs))
end

"""
    dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{T} where { T <: Real }
    dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{Float64} where { T <: Integer }

Compute a spectrum which is `spectrum` rescaled to a live time times probe current equal to `dose`.
Useful for setting spectra on an equivalent acquisition duration scale.
"""
function dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{T} where { T <: AbstractFloat }
    res = copy(spectrum)
    scale = dose / NeXLSpectrum.dose(res)
    res.counts .*= scale
    res[:LiveTime] *= scale
    res[:Name] = "N[$(spectrum[:Name]), $dose nA⋅s]"
    return res
end
function dosenormalize(spectrum::Spectrum{T}, dose=60.0)::Spectrum{Float64} where { T <: Integer }
    scale = dose / NeXLSpectrum.dose(spectrum)
    newProps = copy(spectrum.properties)
    newProps[:LiveTime] *= scale
    newProps[:Name] = "N[$(spectrum[:Name]), $dose nA⋅s]"
    return Spectrum(spectrum.energy, Float64.(spectrum.counts) * scale, newProps)
end

"""
    recalibrate(s::Spectrum{T}, es::LinearEnergyScale)

Allows changing the energy scale on a spectrum from one LinearEnergyScale to another as though the spectrum were 
measured on a different detector.  The algorith uses a FFT-base scheme to rescale and shift the spectral data.
Ths scheme allows for fractional shifts of offset and fractional changes in the width.  It is limited in that
the change in width must produce an integral number of channels in the resulting spectrum.  The algorithm 
maintains the total spectrum integral so the new spectrum can be used for quantitative purposes.
Plotting one spectrum over the other should maintain peak position but is likely to change the channel counts.
"""
function recalibrate(s::Spectrum{T}, es::LinearEnergyScale) where { T <: Real }
    (!(s.energy isa LinearEnergyScale)) && error("The rescale(...) function requires that the spectrum has a LinearEnergyScale.")
    # Resample the spectrum onto a new channel width (offset remains same)
    cxs, oldsum = counts(s), sum(s.counts)
    oldlen, newlen = length(cxs), round(Int, (s.energy.width/es.width)*length(cxs))
    if newlen ≠ length(cxs)
        cxs = FourierTools.resample(cxs, newlen)
    end
    newWidth = s.energy.width*(oldlen/newlen)
    # Shift the first channel from s.energy.offset to es.offset 
    shft = (s.energy.offset - es.offset) / newWidth
    if abs(shft) > 0.01
        cxs = FourierTools.shift(cxs, shft)
    end 
    return Spectrum(es, cxs*(oldsum/sum(cxs)), copy(s.properties))
end

"""
    shift(s::Spectrum, ev::AbstractFloat)::Spectrum

Shift the entire spectrum along the energy axis by a specified number of ev by modifying the
counts data.  This function can shift by a fractional number of channels.
"""
function shift(s::Spectrum, ev::AbstractFloat)::Spectrum
    return Spectrum(s.energy, FourierTools.shift(counts(s), ev/s.energy.width), copy(s.properties))
end

"""
    shannon_entropy(spec::Spectrum)

Computes a measure of the information content in a spectrum.  As there become more and more
distinct values in a spectrum, this value approaches log2(nchannels(spec)).  This number
reflects the number of bits necessary to encode the spectrum data with maximum efficiency.

This is inspired by John Colby's FLAME software which did something similar.  Although, to
be honest, I don't know how his algorithm was implemented.
"""
function shannon_entropy(spec::Spectrum)
    d = Dict{Int,Float64}()
    foreach(c->d[c] = get(d, c, 0.0) + 1.0 / length(spec), counts(spec))
    f(c) = c<=0 ? 0.0 : log2(c)
    return -sum(cx -> cx*f(cx), values(d))
end