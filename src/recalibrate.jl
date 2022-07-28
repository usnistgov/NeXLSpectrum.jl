import FourierTools

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
    (!(s.energy isa LinearEnergyScale)) && error("The recalibrate(...) function requires that the spectrum has a LinearEnergyScale.")
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