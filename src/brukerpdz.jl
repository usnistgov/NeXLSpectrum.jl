using Mmap
using Dates

function readbrukerpdz(fn::AbstractString)
    open(fn,"r") do f
        res = readbrukerpdz(f)
        res[:Name] = splitpath(fn)[end]
        return res
    end
end

"""
    readbrukerpdz(io::IO|fn::AbstractString)

Read the default data file type for the Bruker handheld XRF unit. (Reverse engineered.)
"""
function readbrukerpdz(io::IO)
    mm = Mmap.mmap(io)
    nCh = reinterpret(UInt16, mm[7:8])[1]
    eVpCh = reinterpret(Float64, mm[51:58])[1]
    yr, mon, _, day, hour, min, sec = reinterpret(UInt16,mm[147:160])
    e0, pc = reinterpret(Float32,mm[163:170])
    rt, dt, _, lt = reinterpret(Float32, mm[343:358])
    filt = reinterpret(UInt16,mm[115:130]) # Filter / thickness
    filters=[ (elements[filt[2*i-1]],filt[2*i]*1.0e-6) for i in filter(i->filt[2*i-1]>0, 1:4) ]
    dt = DateTime(yr, mon, day, hour, min, sec)
    props=Dict{Symbol,Any}(:RealTime=>convert(Float64,rt), :LiveTime=>convert(Float64,lt),
            :AcquisitionTime=>dt, :BeamEnergy=>convert(Float64,e0)*1000.0, :ProbeCurrent=>convert(Float64,pc)*1000.0,
            :XRFFilter => filters)
    counts = reinterpret(UInt32, mm[length(mm)-4*nCh+1:length(mm)])
    return Spectrum(LinearEnergyScale(0.0, eVpCh), collect(counts), props)
end
