using Polynomials
using DataFrames
using PeriodicTable

# Keeps track of the number of spectra in this session.

let spectrumIndex = 1
    global spectrumCounter() = (spectrumIndex += 1)
end

"""
    Spectrum
A structure to hold spectrum data (energy scale, counts and metadata). Spectrum implements indexing using various
different mechanisms.  If spec is a Spectrum then

    spec[123] # will return the number of counts in channel 123
	spec[123:223] # will return a Vector of counts from channel 123:223
    spec[134.] # will return the number of counts in the channel at energy 134.0 eV
    spec[:Comment] # will return the property named :Comment

Metadata is identified by a symbol. Predefined symbols include

    :BeamEnergy    # In eV
	:Elevation     # In radians
	:TakeOffAngle  # In radians (Detector position)
    :Azimuthal     # In radians (Detector position)
	:WorkingDistance # In cm
    :LiveTime      # In seconds
    :RealTime      # In seconds
    :ProbeCurrent  # In nano-amps
    :Name          # A string
    :Owner         # A string
    :StagePosition # A Dict{Symbol,Real} with entries :X, :Y, :Z, :R, :T, B: in cm and degrees
    :Comment       # A string
    :Composition   # A Material (known composition, not measured)
    :Elements      # A collection of elements in the material
    :ReferenceROIS # A collection of reference ROIs (as Vector{ReferenceROI})
    :Detector      # A Detector like a SimpleEDS
    :Filename      # Source filename
    :Coating       # A Film (eg. 10 nm of C|Au etc.)
	:AcquisitionTime # Date and time of acquisition (DateTime struct)
	:Signature     # Dict{Element,Real} with the "particle signature"

Less common items

	:ImageMag	   # Magnification (assuming a 3.5" image) of the first image
	:ImageZoom     # Additional zoom for second image in a two image TIFF
	:Operator      # Analyst in ASPEX TIFF files

Not all spectra will define all properties.
"""
struct Spectrum{T<:Real} <: AbstractVector{T}
    energy::EnergyScale
    counts::Vector{T}
    properties::Dict{Symbol,Any}
    hash::UInt  # random number (stays fixed as underlying data changes)

    function Spectrum(energy::EnergyScale, data::Vector{<:Real}, props::Dict{Symbol,Any})
        props[:Name] = get(props, :Name, "Spectrum[$(spectrumCounter())]")
        return new{typeof(data[1])}(energy, data, props, _hashsp(energy, data, props))
    end
end

_hashsp(e, d, p) = xor(hash(e), xor(hash(d), hash(p)))

Base.hash(spec::Spectrum, h::UInt) = hash(spec.hash, h)
Base.isequal(spec1::Spectrum, spec2::Spectrum) =
    (hash(spec1) == hash(spec2)) &&
    isequal(spec1.energy, spec2.energy) &&
    isequal(spec1.properies, spec2.properties) && isequal(spec1.counts, spec2.counts)
Base.isless(s1::Spectrum, s2::Spectrum) =
    isequal(s1[:Name], s2[:Name]) ? isless(s1.hash, s2.hash) : isless(s1[:Name], s2[:Name])
# Make it act like an AbstractVector
Base.eltype(spec::Spectrum) = Spectrum
Base.length(spec::Spectrum) = length(spec.counts)
Base.ndims(spec::Spectrum) = ndims(spec.counts)
Base.size(spec::Spectrum) = size(spec.counts)
Base.size(spec::Spectrum, n) = size(spec.counts, n)
Base.axes(spec::Spectrum) = axes(spec.counts)
Base.axes(spec::Spectrum, n) = axes(spec.counts, n)
Base.eachindex(spec::Spectrum) = eachindex(spec.counts)
Base.stride(spec::Spectrum, k) = stride(spec.counts, k)
Base.strides(spec::Spectrum) = strides(spec.counts)
Base.getindex(spec::Spectrum, idx::Int) = spec.counts[idx]
Base.getindex(spec::Spectrum, sr::StepRange{Int64,Int64}) = spec.counts[sr]
Base.getindex(spec::Spectrum, ur::UnitRange{Int64}) = spec.counts[ur]
Base.get(spec::Spectrum, idx::Int, def = convert(typeof(spec.counts[1]), 0)) = get(spec.counts, idx, def)
Base.setindex!(spec::Spectrum, val::Real, idx::Int) = spec.counts[idx] = convert(typeof(spec.counts[1]), val)
Base.setindex!(spec::Spectrum, vals, ur::UnitRange{Int}) = spec.counts[ur] = vals
Base.setindex!(spec::Spectrum, vals, sr::StepRange{Int}) = spec.counts[sr] = vals
Base.copy(spec::Spectrum) = Spectrum(spec.energy, copy(spec.counts), copy(spec.properties))

rangeofenergies(spec::Spectrum, ch) = (energy(ch, spec.energy), energy(ch + 1, spec.energy))

"""
	minrequired(::Type{XXX})

Returns the minimum required properties.  Other classes implement this to check whether a Spectrum has all the
necessary properties for the specified algorithm or data structure.
"""
minproperties(::Type{Spectrum}) = (:Name)


"""
    hasminrequired(ty::Type, spec::Spectrum)

Does this spectrum have the minimal set of required properties?
"""
hasminrequired(ty::Type, spec::Spectrum) = all(haskey(spec.properties, a) for a in minproperties(ty))

"""
    requiredbutmissing(ty::Type, spec::Spectrum)

List any required but missing properties.
"""
requiredbutmissing(ty::Type, spec::Spectrum) = filter(a -> !haskey(spec.property, a), minproperties(ty))


maxproperty(specs, prop::Symbol) = maximum(spec -> spec[prop], specs)
minproperty(specs, prop::Symbol) = minimum(spec -> spec[prop], specs)
sameproperty(specs, prop::Symbol) = all(spec -> spec[prop] == specs[1][prop], specs) ? #
    specs[1][prop] : #
    error("The property $prop is not equal for all these spectra.")

function Base.show(io::IO, spec::Spectrum)
    comp = haskey(spec, :Composition) ? name(spec[:Composition]) : "Unknown"
    e0 = haskey(spec, :BeamEnergy) ? "$(round(spec[:BeamEnergy]/1000.0,sigdigits=3)) keV" : "Unknown keV"
    print(io, "Spectrum[$(spec[:Name]), $e0, $comp, $(round(sum(spec.counts),sigdigits=3)) counts]")
end

function textplot(io::IO, spec::Spectrum; size = (16, 80))
    (rows, cols) = size
    # how much to plot
    e0_eV = haskey(spec, :BeamEnergy) ? spec[:BeamEnergy] : energy(length(spec), spec)
    maxCh = min(channel(e0_eV, spec), length(spec))
    step, max = maxCh ÷ cols, maximum(spec.counts)
    maxes = [rows * (maximum(spec.counts[(i-1)*step+1:i*step]) / max) for i = 1:cols]
    for r = 1:rows
        ss = ""
        for i = 1:cols
            ss = ss * (r ≥ rows - maxes[i] ? "*" : " ")
        end
        if r == 1
            println(io, "$ss $max")
        elseif r == rows
            println(io, "$ss $(0.001*e0_eV) keV]")
        else
            println(io, ss)
        end
    end
end

"""
    asa(::Type{DataFrame}, spec::Spectrum)

Converts the spectrum energy and counts data into a DataFrame.
"""
NeXLUncertainties.asa(::Type{DataFrame}, spec::Spectrum; properties::Bool = false)::DataFrame =
    properties ? DataFrame(Keys = [keys(spec.properties)...], Values = repr.([values(spec.properties)...])) : #
    DataFrame(E = energyscale(spec), I = counts(spec))

"""
    apply(spec::Spectrum, det::SimpleEDS)::Spectrum

Applys the specified detector to this spectrum by ensuring that the energy scales match and spec[:Detector] = det.
Creates a copy of the original spectrum unless the detector is already det.
"""

function apply(spec::Spectrum, det::SimpleEDS)::Spectrum
    if (haskey(spec, :Detector)) && (spec[:Detector] == det)
        return spec
    else
        props = copy(spec.properties)
        props[:Detector] = det
        return Spectrum(det.scale, spec.counts, props)
    end
end

"""
	matching(spec::Spectrum, res::Resolution, lld::Int=1)::EDSDetector

Build an EDSDetector to match the channel count and energy scale in this spectrum.
"""
matching(spec::Spectrum, res::Resolution, lld::Int = 1)::EDSDetector = SimpleEDS(spec, res, lld)

"""
	matching(spec::Spectrum, resMnKa::Float64, lld::Int=1)::SimpleEDS

Build an EDSDetector to match the channel count and energy scale in this spectrum.
"""
matching(spec::Spectrum, resMnKa::Float64, lld::Int = 1)::SimpleEDS =
    SimpleEDS(length(spec), spec.energy, MnKaResolution(resMnKa), lld)

matching(spec::Spectrum, resMnKa::Float64, lld::Int, minByFam::Dict{Shell,Element})::BasicEDS =
    BasicEDS(length(spec), spec.energy, MnKaResolution(resMnKa), lld, minByFam)

setproperty!(spec::Spectrum, sym::Symbol, val::Any) = setindex!(props, sym, val)
Base.get(spec::Spectrum, sym::Symbol, def::Any = missing) = get(spec.properties, sym, def)
Base.getindex(spec::Spectrum, sym::Symbol)::Any = spec.properties[sym]
Base.setindex!(spec::Spectrum, val::Any, sym::Symbol) = spec.properties[sym] = val

Base.keys(spec::Spectrum) = keys(spec.properties)

Base.haskey(spec::Spectrum, sym::Symbol) = haskey(spec.properties, sym)

"""
    elms(spec::Spectrum, withcoating = false, def=missing)

Returns a list of the elements associated with this spectrum. `withcoating` determines whether the coating
elements are also added.
"""
function NeXLCore.elms(spec::Spectrum, withcoating = false, def = missing)
    res = Set{Element}()
    if haskey(spec, :Elements)
        append!(spec[:Elements])
    elseif haskey(spec, :Composition)
        append!(keys(spec[:Composition]))
    elseif haskey(spec, :Signature)
        append!(keys(spec[:Signature]))
    end
    if withcoating && haskey(spec, :Coating)
        append!(res, keys(material(spec[:Coating])))
    end
    return length(res) == 0 ? def : res
end

"""
    dose(spec::Spectrum, def=missing)

The probe dose in nano-amp seconds
"""
function dose(spec::Spectrum, def = missing)::Union{Float64,Missing}
    res = get(spec.properties, :LiveTime, missing) * get(spec, :ProbeCurrent, missing)
    return isequal(res, missing) ? def : res
end

"""
    NeXLCore.energy(ch::Int, spec::Spectrum)

The energy of the start of the ch-th channel.
"""
NeXLCore.energy(ch::Int, spec::Spectrum)::Float64 = energy(ch, spec.energy)

"""
    channelwidth(ch::Int, spec::Spectrum)::Float64

Returns the width of the `ch` channel
"""
channelwidth(ch::Int, spec::Spectrum) = energy(ch + 1, spec) - energy(ch, spec)

"""
    channel(eV::Float64, spec::Spectrum)

The index of the channel containing the specified energy.
"""
channel(eV::Float64, spec::Spectrum)::Int = channel(eV, spec.energy)

"""
    counts(spec::Spectrum, numType::Type{T}, applyLLD=false)::Vector{T} where {T<:Number}

Creates a copy of the spectrum counts data as the specified Number type. If the spectrum has a :Detector
property then the detector's lld (low-level discriminator) and applyLLD=true then the lld is applied to the result
by setting all channels less-than-or-equal to det.lld to zero.
"""
function counts(spec::Spectrum, numType::Type{T} = Float64, applyLLD = false) where {T<:Number}
    res = map(n -> convert(numType, n), spec.counts)
    if applyLLD && haskey(spec, :Detector)
        applylld(spec, lld(spec[:Detector]))
    end
    return res
end

function applylld(spec::Spectrum, lld::Int)
    fill!(view(spec, 1:lld), spec[lld+1])
end

"""
    counts(spec::Spectrum, channels::UnitRange{Int}, numType::Type{T}, applyLLD=false)::Vector{T} where {T<:Number}

Creates a copy of the spectrum counts data as the specified Number type.  If the spectrum has a :Detector
property then the detector's lld (low-level discriminator) and applyLLD=true then the lld is applied to the result
by setting all channels less-than-or-equal to det.lld to zero.
"""
function counts(spec::Spectrum, channels::UnitRange{Int}, numType::Type{T}, applyLLD = false) where {T<:Real}
    res = map(n -> convert(numType, n), spec.counts[channels])
    if applyLLD && haskey(spec, :Detector)
        fill!(view(res, 1:lld(spec)-channels.start+1), zero(numType))
    end
    return res
end

"""
    lld(spec::Spectrum)

Gets the low-level discriminator associated with this spectrum if there is one.
"""
lld(spec::Spectrum) = haskey(spec.properties, :Detector) ? lld(spec.properties[:Detector]) : 1


"""
    normalizeDoseWidth(spec::Spectrum)::Vector{Float64}

Normalize the channel intensities to counts/(nA⋅s⋅eV).  Good for comparing spectra collected at different detector
channel widths.
"""
function normalizedosewidth(spec::Spectrum, defDose = missing)::Vector{Float64}
    ds = dose(spec, defDose)
    if ismissing(ds)
        error("The required spectrum dose in not available in normalizeDoseWidth(spec).")
    end
    return map(ch -> convert(Float64, spec.counts[ch]) / (ds * channelwidth(ch, spec)), eachindex(spec.counts))
end

"""
	findmax(spec::Spectrum, chs::UnitRange{Int})

Returns the (maximum intensity, channel index) over the specified range of channels
"""
function Base.findmax(spec::Spectrum, chs::UnitRange{Int})
    max = findmax(spec.counts[chs])
    return (max[1] + chs.start - 1, max[2])
end

"""
	findmax(spec::Spectrum)

Returns the (maximum intensity, channel index) over all channels
"""
Base.findmax(spec::Spectrum) = findmax(spec.counts)

"""
   integrate(spec, channels)

Sums all the counts in the specified channels.  No background correction.
"""
integrate(spec::Spectrum, channels::UnitRange{Int}) = sum(spec.counts[channels])

"""
   integrate(spec, channels)

Sums all the counts in the specified energy range.  No background correction.
"""
integrate(spec::Spectrum, energyRange::StepRangeLen{Float64}) =
    sum(spec.counts[channel(energyRange[1], spec):channel(energyRange[end], spec)])

"""
    integrate(spec::Spectrum, back1::UnitRange{Int}, peak::UnitRange{Int}, back2::UnitRange{Int})::Float64
Perform a background corrected peak integration using channel ranges. Fits a line to each background region and
extrapolates the background from the closest background channel through the peak region.
"""
function integrate(spec::Spectrum, back1::UnitRange{Int}, peak::UnitRange{Int}, back2::UnitRange{Int})::Float64
    p1 = Polynomials.fit(Poly, back1, spec.counts[back1], 1)
    p2 = Polynomials.fit(Poly, back2, spec.counts[back2], 1)
    c1, c2 = back1.stop, back2.start
    i1, i2 = p1(c1), p2(c2)
    m = (i2 - i1) / (c2 - c1)
    back = Poly([i1 - m * c1, m])
    sum(spec.counts[peak] - map(back, peak))
end

"""
    integrate(spec::Spectrum, back1::StepRangeLen{Float64}, peak::StepRangeLen{Float64}, back2::StepRangeLen{Float64})::Float64
Perform a background corrected peak integration using energy (eV) ranges. Converts the energy ranges to channels
ranges before performing the integral.
"""
function integrate(
    spec::Spectrum,
    back1::StepRangeLen{Float64},
    peak::StepRangeLen{Float64},
    back2::StepRangeLen{Float64},
)::Float64
    b1i = channel(back1[back1.offset], spec):channel(back1[end], spec)
    b2i = channel(back2[back2.offset], spec):channel(back2[end], spec)
    pi = channel(peak[peak.offset], spec):channel(peak[end], spec)
    integrate(spec, b1i, pi, b2i)
end

"""
    integrate(spec::Spectrum)

Total integral of all counts from the LLD to the beam energy
"""
function integrate(spec::Spectrum)
    last = min(haskey(spec, :BeamEnergy) ? channel(spec[:BeamEnergy], spec) : length(spec.counts), length(spec.counts))
    return integrate(spec, lld(spec):last)
end

"""
    findmax(spec::Spectrum)
"""
function findmax(spec::Spectrum)
    last = min(haskey(spec, :BeamEnergy) ? channel(spec[:BeamEnergy], spec) : length(spec.counts), length(spec.counts))
    return findmax(spec.counts, lld(spec):last)
end

"""
    energyscale(spec::Spectrum)

Returns an array with the bin-by-bin energies
"""
energyscale(spec::Spectrum) = energyscale(spec.energy, 1:length(spec))

"""
    simpleEDS(spec::Spectrum, fwhmatmnka::Float64)

Build a SimpleEDS object for this spectrum with the specified FWHM at Mn Kα.
"""
simpleEDS(spec::Spectrum, fwhmatmnka::Float64) = SimpleEDS(length(spec), spec.energy, MnKaResolution(fwhmatmnka))


"""
    subsample(spec::Spectrum, frac::Float64)

Subsample the counts data in a spectrum according to a statistically valid algorithm.  Returns
`spec` if frac>=1.0.
"""
function subsample(spec::Spectrum, frac::Float64)::Spectrum
    @assert frac > 0.0 "frac must be larger than zero."
    @assert frac <= 1.0 "frac must be less than or equal to 1.0."
    @assert haskey(spec, :LiveTime) "Please specify a :LiveTime in subsample"
    ss(n, f) = n > 0 ? sum(rand() <= f ? 1 : 0 for _ = 1:n) : 0
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
function subdivide(spec::Spectrum, n::Int)::Vector{Spectrum}
    res = zeros(Int, n, length(spec.counts))
    for ch in eachindex(spec.counts)
        # Assign each event to one and only one detector
        si = rand(1:n, floor(Int, spec[ch]))
        for i = 1:n
            res[i, ch] = count(e -> e == i, si)
        end
    end
    specs = Spectrum[]
    for i = 1:n
        props = deepcopy(spec.properties)
        props[:Name] = "Sub[$(spec[:Name]),$(i) of $(n)]"
        props[:LiveTime] = spec[:LiveTime] / n # Must have
        if haskey(spec, :RealTime)
            props[:RealTime] = spec[:RealTime] / n # Might have
        end
        push!(specs, Spectrum(spec.energy, res[i, :], props))
    end
    return specs
end


"""
    estimatebackground(data::AbstractArray{Float64}, channel::Int, width::Int=5, order::Int=2)

Returns the tangent to the a quadratic fit to the counts data centered at channel with width
"""
function estimatebackground(data::AbstractArray{Float64}, channel::Int, width::Int = 5, order::Int = 2)::Poly
    minCh, maxCh = max(1, channel - width), min(length(data), channel + width)
    if maxCh - minCh >= order
        fit = Polynomials.fit(Poly, (minCh-channel):(maxCh-channel), data[minCh:maxCh], order)
        return Poly([fit(0), polyder(fit)(0)]) # Linear
    else
        return Poly([mean(data[minCh:maxCh]), 0.0])
    end
end

"""
    modelBackground(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)

spec: A spectrum containing a peak centered on chs
chs:  A range of channels containing a peak
ash:  The edge (as an AtomicSubShell)

A simple model for modeling the background under a characteristic x-ray peak. The model
fits a line to low and high energy background regions around chs.start and chs.end. If
the low energy line extended out to the edge energy is larger than the high energy line
at the same energy, then a negative going edge is fit between the two. Otherwise a line
is fit between the low energy side and the high energy side. This model only works when
there are no peak interference over the range chs.
"""
function modelBackground(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)
    cnts, ec = counts(spec), channel(energy(ash), spec)
    bl, bh = estimatebackground(cnts, chs.start, 5), estimatebackground(cnts, chs.stop, 5)
    # bh = ch-> mean(cnts[chs.stop:min(length(cnts),chs.stop+5)])
    if (ec < chs.stop) && (bl(ec - chs.start) > bh(ec - chs.stop)) && (energy(ash) < 2.0e3)
        res = zeros(Float64, length(chs))
        for y = chs.start:ec-1
            res[y-chs.start+1] = bl(y - chs.start)
        end
        res[ec-chs.start+1] = 0.5 * (bl(ec - chs.start) + bh(ec - chs.stop))
        for y = ec+1:chs.stop
            res[y-chs.start+1] = bh(y - chs.stop)
        end
    else
        s = (bh(0) - bl(0)) / length(chs)
        back = Poly([bl(0), s])
        res = back.(collect(0:length(chs)-1))
    end
    return res
end


"""
    modelBackground(spec::Spectrum, chs::UnitRange{Int})

spec: A spectrum containing a peak centered on chs
chs:  A range of channels containing a peak

A simple model for modeling the background under a characteristic x-ray peak. The model
fits a line between the  low and high energy background regions around chs.start and chs.end.
This model only works when there are no peak interference over the range chs.
"""
function modelBackground(spec::Spectrum, chs::UnitRange{Int})
    bl = estimatebackground(counts(spec), chs.start, 5)
    bh = estimatebackground(counts(spec), chs.stop, 5)
    s = (bh(0) - bl(0)) / length(chs)
    back = Poly([bl(0), s])
    return back.(collect(0:length(chs)-1))
end

"""
    extractcharacteristic(spec::Spectrum, lowBack::UnitRange{Int}, highBack::UnitRange{Int})::Vector{Float64}

Extract the characteristic intensity for the peak located within chs with an edge at ash.
"""
function extractcharacteristic(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)::Vector{Float64}
    return counts(spec, chs, Float64) - modelBackground(spec, chs, ash)
end


"""
    peak(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)::Float64

Estimates the peak intensity for the characteristic X-ray in the specified range of channels.
"""
peak(spec::Spectrum, chs::UnitRange{Int})::Float64 = return sum(counts(spec, chs, Float64)) - back(spec, chs)

"""
    back(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)::Float64

Estimates the background intensity for the characteristic X-ray in the specified range of channels.
"""
back(spec::Spectrum, chs::UnitRange{Int})::Float64 = sum(modelBackground(spec, chs))


"""
    peaktobackground(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)::Float64

Estimates the peak-to-background ratio for the characteristic X-ray intensity in the specified range of channels
which encompass the specified AtomicSubShell.
"""
function peaktobackground(spec::Spectrum, chs::UnitRange{Int}, ash::AtomicSubShell)::Float64
    back = sum(modelBackground(spec, chs, ash))
    return (sum(counts(spec, chs, Float64)) - back) / back
end

"""
    estkratio(unk::Spectrum, std::Spectrum, chs::UnitRange{Int})

Estimates the k-ratio from niave models of peak and background intensity.  Only works if the peak is not interfered.
"""
estkratio(unk::Spectrum, std::Spectrum, chs::UnitRange{Int}) = peak(unk, chs) * dose(std) / (peak(std, chs) * dose(unk))


"""
    NeXLUncertainties.asa(::Type{DataFrame}, spec::AbstractVector{Spectrum})::DataFrame

Returns a DataFrame that summarizes the list of spectra.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, specs::AbstractVector{Spectrum})::DataFrame
    _asname(comp) = ismissing(comp) ? missing : name(comp)
    unf, unl, uns = Union{Float64,Missing}, Union{Film,Nothing}, Union{String,Missing}
    nme, e0, pc, lt, rt, coat, integ, comp = String[], unf[], unf[], unf[], unf[], unl[], Float64[], uns[]
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
function NeXLUncertainties.asa(::Type{DataFrame}, stds::AbstractDict{Element,Spectrum})::DataFrame
    _asname(comp) = ismissing(comp) ? missing : name(comp)
    unf, unl, uns = Union{Float64,Missing}, Union{Film,Nothing}, Union{String,Missing}
    elm, zs, mfs, nme, e0, pc, lt, rt, coat, integ, comp =
        String[], Int[], Float64[], String[], unf[], unf[], unf[], unf[], unl[], Float64[], uns[]
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

Outputs a description of the data in the spectrum.
"""
function details(io::IO, spec::Spectrum)
    println(io, "           Name:   $(spec[:Name])")
    println(io, "    Beam energy:   $(get(spec, :BeamEnergy, missing)/1000.0) keV")
    println(io, "  Probe current:   $(get(spec, :ProbeCurrent, missing)) nA")
    println(io, "      Live time:   $(get(spec, :LiveTime, missing)) s")
    println(io, "        Coating:   $(get(spec,:Coating, "None"))")
    println(io, "       Detector:   $(get(spec, :Detector, missing))")
    println(io, "        Comment:   $(get(spec, :Comment, missing))")
    println(io, "       Integral:   $(integrate(spec)) counts")
    comp = get(spec, :Composition, missing)
    if !ismissing(comp)
        println(io, "    Composition:   $(comp)")
        det = get(spec, :Detector, missing)
        if !ismissing(det)
            coating = get(spec, :Coating, missing)
            comp2 = collect(keys(comp))
            if !ismissing(coating)
                append!(comp2, keys(coating.material))
            end
            for elm1 in keys(comp)
                for ext1 in extents(elm1, det, 1.0e-4)
                    print(io, "            ROI:   $(elm1.symbol)[$(ext1)]")
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
                        println(io, " intersects $(join(intersects,", "))")
                    else
                        p, b = peak(spec, ext1), back(spec, ext1)
                        σ = p / sqrt(b)
                        println(io, " = $(round(Int,p)) counts over $(round(Int,b)) counts - σ = $(round(Int,σ))")
                    end
                end
            end
        end
    end
    return nothing
end

"""
    details(spec::Spectrum)

Outputs a description of the data in the spectrum to standard output.
"""
details(spec::Spectrum) = details(stdout, spec)

function commonproperties(props1::Dict{Symbol,Any}, props2::Dict{Symbol,Any})
    res = Dict{Symbol,Any}()
    for (key, sp1) in props1
        if isequal(sp1, get(props2, key, missing))
            res[key] = sp1
        end
    end
    return res
end

commonproperties(specs::AbstractArray{Spectrum}) =
    reduce((props, sp2) -> commonproperties(props, sp2.properties), specs[2:end], specs[1].properties)

function maxspectrum(spec1::Spectrum, spec2::Spectrum)
    maxi(a, b) = collect(max(a[i], b[i]) for i in eachindex(a))
    props = commonproperties(spec1.properties, spec2.properties)
    props[:Name] = "MaxSpectrum"
    return Spectrum(spec1.energy, maxi(counts(spec1), counts(spec2)), props)
end

maxspectrum(specs::AbstractArray{Spectrum}) = reduce(maxspectrum, specs)

maxspectrum(specs::Vararg{Spectrum}) = reduce(maxspectrum, specs)
