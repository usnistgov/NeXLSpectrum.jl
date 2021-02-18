"""
SpectrumScaling types are designed to rescale spectrum data primarily for plotting.

**Implement**

    Base.show(io::IO, scn::SpectrumScaling)
    scalefactor(sc::SpectrumScaling, spec::Spectrum)
"""
abstract type SpectrumScaling end

scaledcounts(sc::SpectrumScaling, spec::Spectrum) = scalefactor(sc, spec) * counts(spec)

"""
Don't scale the spectrum data.
"""
struct NoScaling <: SpectrumScaling end
Base.show(io::IO, scn::NoScaling) = print(io, "Counts")
scalefactor(sc::NoScaling, spec::Spectrum) = 1.0

"""
Scale to a constant dose⋅width (Counts/(nA⋅s/eV))  Requires a spectrum has both :ProbeCurrent & :LiveTime.
"""
struct ScaleDoseWidth <: SpectrumScaling
    defDose::Float64  # default in nA⋅s/eV

    ScaleDoseWidth(nAspeV::Float64 = 60.0) = new(nAspeV)
end

Base.show(io::IO, sc::ScaleDoseWidth) = print(io, "Counts/($(sc.defDose) nA⋅s/eV)")

scalefactor(sc::ScaleDoseWidth, spec::Spectrum) =
    sc.defDose / (dose(spec) * channelwidth(1, spec))

function scaledcounts(sc::ScaleDoseWidth, spec::Spectrum)
    # Specialized since each channel can have a different width and thus scale...
    ds = dose(spec)
    return map(ch -> spec[ch] * sc.defDose / (ds * channelwidth(ch, spec)), eachindex(spec))
end

"""
Scale to a constant dose⋅width (Counts/(nA⋅s/eV))  Requires a spectrum has both :ProbeCurrent & :LiveTime.
"""
struct ScaleWidth <: SpectrumScaling end

Base.show(io::IO, sc::ScaleWidth) = print(io, "Counts/eV")

scalefactor(sc::ScaleWidth, spec::Spectrum) = channelwidth(1, spec)

# Specialized since each channel can have a different width and thus scale...
scaledcounts(sc::ScaleWidth, spec::Spectrum) =
    map(ch -> spec[ch] / channelwidth(ch, spec), eachindex(spec))


"""
Scale to a constant dose (Counts/(nA⋅s)).   Requires a spectrum has both :ProbeCurrent & :LiveTime.
"""
struct ScaleDose <: SpectrumScaling
    defDose::Float64  # default in nA⋅s

    ScaleDose(nAs::Float64 = 60.0) = new(nAs)
end

Base.show(io::IO, sc::ScaleDose) = print(io, "Counts/($(sc.defDose) nA⋅s)")
scalefactor(sc::ScaleDose, spec::Spectrum) = sc.defDose / dose(spec)

"""
Scale to a fixed total integral.
"""
struct ScaleSum <: SpectrumScaling
    defSum::Float64  # default in counts

    ScaleSum(counts::Float64 = 1.0e5) = new(counts)
end

Base.show(io::IO, sc::ScaleSum) = print(io, "Total Integral ($(sc.defSum))")
scalefactor(sc::ScaleSum, spec::Spectrum) = sc.defSum / integrate(spec)


"""
Scale to a fixed peak intensity
"""
struct ScalePeak <: SpectrumScaling
    defPeak::Float64  # default in counts

    ScalePeak(peak::Float64 = 1000.0) = new(peak)
end
Base.show(io::IO, sc::ScalePeak) = print(io, "Peak Counts ($(sc.defPeak))")
scalefactor(sc::ScalePeak, spec::Spectrum) =
    sc.defPeak / Base.findmax(spec[lld(spec):end])[1]


"""
Scale to a default sum in the specified ROI.
"""
struct ScaleROISum <: SpectrumScaling
    roi::UnitRange{Int}
    defSum::Float64

    ScaleROISum(roi::UnitRange{Int}, defSum::Float64) = new(roi)
end

Base.show(io::IO, sc::ScaleROISum) = print(io, "ROI Integral [$(sc.roi), $(sc.defSum)]")
scalefactor(sc::ScaleROISum, spec::Spectrum) = sc.defSum / findmax(sp[roi])[1]
