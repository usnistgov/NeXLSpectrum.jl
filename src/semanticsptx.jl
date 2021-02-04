using OnlineStats: Mean, value, fit!

"""
    readptx(fn::AbstractString, scale::EnergyScale, nch::Int, dets=[true,true,true,true],frames=1:1000)

Read a SEMantics .ptx file into a HyperSpectrum.  Really needs a better ZIP file reader to work....
"""
function readptx(
    fn::AbstractString,
    scale::EnergyScale,
    nch::Int,
    dets = [true, true, true, true],
    frames = 1:1000,
)::HyperSpectrum
    function readmsg(io::IO)
        hdr = Array{UInt8}(undef, 32)
        if readbytes!(io, hdr) == 32
            name = String(hdr[1:findfirst(c -> c == 0x0, hdr)-1])
            len = reinterpret(UInt32, hdr[17:20])[1]
            msg = Array{UInt8}(undef, len)
            if readbytes!(io, msg) == len
                return (name, msg)
            end
        end
        return ("End", nothing)
    end
    zip = ZipFile.Reader(fn)
    hdrfile = zip.files[findfirst(f -> f.name == "Header", zip.files)]
    xml = read(hdrfile, String)
    hdr = readxml(IOBuffer(xml))
    props = Dict{Symbol,Any}()
    dwell = 1.0e-6 * parse(Float64, findfirst("/OpticState/mDwell", hdr).content)
    fov = 0.1 * parse(Float64, findfirst("/OpticState/mFieldOfView", hdr).content)
    width = parse(Int, findfirst("/OpticState/mWidth", hdr).content)
    height = parse(Int, findfirst("/OpticState/mHeight", hdr).content)
    rx = parse(Int, findfirst("/OpticState/mRaster/x", hdr).content)
    ry = parse(Int, findfirst("/OpticState/mRaster/y", hdr).content)
    rw = parse(Int, findfirst("/OpticState/mRaster/width", hdr).content)
    rh = parse(Int, findfirst("/OpticState/mRaster/height", hdr).content)
    # The output array consists of list of datagrams, each has got five double-word (32-bit) values (DW0 to
    # DW4). Little endian is storage is applied. Most significant byte (bits 31-24) of W0 tells the event type:
    # W0 = 0x00xxxxxx – measured energy
    # W0 – bits 0-7: filter index (1= filter 1, 2 = filter 2, …)  bits 8-11: detector index (0=1st detector, 1=2nd detector, ...)
    # W1 – energy in eV (24.8 signed fixed point number)
    # W2 – timestamp in 10ns
    # W3 – pixel index
    # W4 – frame index
    # W0 = 0x10xxxxxx – line sync
    # W0 – bit 0: “Line” signal state bits 8-11: detector index (0=1st detector, 1=2nd detector, ...)
    # W2 – timestamp in 10ns
    # W4 – frame index
    # W0 = 0x11xxxxxx – pileup
    # W0 – bits 8-11: detector index (0=1st detector, 1=2nd detector, …)
    # W2 – timestamp in 10ns
    # W4 – frame index
    # W0 = 0x13xxxxxx – preamplifier reset
    # W0 – bits 8-11: detector index (0=1st detector, 1=2nd detector, …)
    # W2 – timestamp in 10ns
    # W4 – frame index
    # W0 = 0x80xxxxxx – dead-time measurement
    # W0 – bits 0-7: filter index (1= filter 1, 2 = filter 2, …) bits 8-11: detector index (0=1st detector, 1=2nd detector, …)
    # W1 – real-time timer value in 10 ns
    # W2 – dead-time timer value in 10 ns
    # W3 – event counter
    hss = HyperSpectrum(scale, props, (rh - ry, rw - rx), nch, UInt16)
    data = zip.files[findfirst(f -> f.name == "Data", zip.files)]
    real, dead, deadpct = zeros(UInt64, length(dets)), zeros(UInt64, length(dets)), Mean()
    seen = Set{String}()
    while true
        cmd, msg = readmsg(data)
        # println("$cmd[$(isnothing(msg) ? "nothing" : length(msg))]")
        if cmd == "End"
            break
        elseif cmd == "GetMapData"
            err = reinterpret(Int32, msg[1:4])
            overflow = reinterpret(UInt32, msg[5:8])
            usage = reinterpret(Float32, msg[9:12])
            encoding = reinterpret(UInt32, msg[13:16])
            ws = reinterpret(UInt32, msg[13:end])
            for i = 1:length(ws)÷5
                w0, w1, w2, w3, w4 = ws[5*(i-1)+1:5*i]
                # println("$(ws[5*(i-1)+1:5*i])")
                x_det = 1 + ((w0 >>> 8) & 0x0F)
                # msg[2] # timestamp in 10 ns units
                msgid = (w0 >>> 48) & 0xFF
                if msgid == 0x80 # dead time
                    if x_det <= length(dets) && dets[x_det] && w4 in frames
                        # println("real=>$w1, dead=>$w2")
                        real[x_det] += w1
                        dead[x_det] += w2
                        fit!(deadpct, convert(Float64, w2) / convert(Float64, w1))
                    end
                elseif msgid == 0x00 # measured energy
                    if x_det <= length(dets) && dets[x_det] && w4 in frames
                        ch = channel(convert(Float64, w1 ÷ 256), scale)
                        r, c = w3 ÷ rw - ry + 1, w3 % rw - rx + 1
                        # println("r,c=>($r, $c) => $ch @ $(w1÷256)")
                        if ch >= 1 && ch <= nch && r >= 1 && r <= rw && c >= 1 && c <= rh
                            hss[r, c, ch] += 1
                        end
                    end
                end
            end
        else
            if !(cmd in seen)
                println("Skipping: $cmd")
                push!(seen, cmd)
            end
        end
    end
    hss[:Dwell] = dwell
    hss[:Elapse] = 1.0e-8 * real # second
    hss[:RealTime] = 1.0e-8 * real / ((rh - ry) * (rw - rx))
    hss[:LiveTime] = 1.0e-8 * (real - dead) / ((rh - ry) * (rw - rx)) # seconds
    hss[:DeadPercent] = OnlineStats.value(deadpct)
    return hss
end
