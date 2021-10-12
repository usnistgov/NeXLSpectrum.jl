using NeXLSpectrum

path=raw"C:\Users\nritchie\OneDrive - National Institute of Standards and Technology (NIST)\Documents\FY2016\SE1610x\Manual\Particle 12 keV"
p = joinpath(path,"map.ptx")

props = Dict(:BeamEnergy => 12.0e3, :ProbeCurrent=>1.0, :RealTime => 0.129202)

ht1 = @time readptx(p, LinearEnergyScale(0.0,10.0), props, 1024, blocksize=8, dets=(true, false, false, false))
ht2 = @time readptx(p, LinearEnergyScale(0.0,10.0), props, 1024, blocksize=8, dets=(false, true, false, false))
ht3 = @time readptx(p, LinearEnergyScale(0.0,10.0), props, 1024, blocksize=8, dets=(false, false, true, false))
ht4 = @time readptx(p, LinearEnergyScale(0.0,10.0), props, 1024, blocksize=8, dets=(false, false, false, true))

specs = [ ht1, ht2, ht3, ht4 ]

NeXLSpectrum.checkmultispec(specs)

map(h->h[:RealTime], specs)