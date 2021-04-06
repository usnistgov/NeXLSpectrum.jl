using p7zip_jll
using Colors

"""
    readptx(
        fn::AbstractString, 
        scale::EnergyScale, # Energy scale for hyperspectrum
        nch::Int, # Number of channels in hyperspectrum
        blocksize = 1,  # Size of blocks to sum to create averaged hyperspectrum
        dets=[true,true,true,true], # Which detectors to include
        frames=1:typemax(Int)) # which frames to include in hyperspectrum

Read a SEMantics .ptx file into a HyperSpectrum.
"""
function readptx(
    fn::AbstractString,
    scale::EnergyScale,
    nch::Int,
    blocksize = 1,
    frames = 1:typemax(Int),
    dets = [true, true, true, true],
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
    mktempdir() do tmp_dir
        p7zip() do p7z_exe
            # Use pkzip to unzip the contents to a temporary directory
            run(`$p7z_exe -bsp0 -bso0 e "$fn" -o"$tmp_dir"`)
        end
        # Read and interpret the header
        xml = read(joinpath(tmp_dir,"Header")) 
        hdr = readxml(IOBuffer(xml))
        props = Dict{Symbol,Any}()
        props[:Name] = splitpath(fn)[end]
        dwell = 1.0e-6 * parse(Float64, findfirst("/OpticState/mDwell", hdr).content)
        fov = 0.1 * parse(Float64, findfirst("/OpticState/mFieldOfView", hdr).content)
        width = parse(Int, findfirst("/OpticState/mWidth", hdr).content)
        height = parse(Int, findfirst("/OpticState/mHeight", hdr).content)
        rx = parse(Int, findfirst("/OpticState/mRaster/x", hdr).content)
        ry = parse(Int, findfirst("/OpticState/mRaster/y", hdr).content)
        rw = parse(Int, findfirst("/OpticState/mRaster/width", hdr).content)
        rh = parse(Int, findfirst("/OpticState/mRaster/height", hdr).content)
        hss = HyperSpectrum(scale, props, ((rh - ry)÷blocksize, (rw - rx)÷blocksize), nch, UInt32)

        open(joinpath(tmp_dir, "Data")) do data
            realtime, deadtime = zero(UInt64), zero(UInt64)
            maxchannel, maxframeid, minframeid = zero(Int32), zero(Int32), typemax(Int32)
            # elapsetime, prevelapse = zero(UInt64), typemax(UInt32)
            seen = Set{String}()
            images = Dict{Tuple{Int32,Int32},Array}()
            cmd, msg = readmsg(data)
            while cmd != "End"
                # println("$cmd[$(isnothing(msg) ? "nothing" : length(msg))]")
                if cmd == "GetMapData"
                    # The output array consists of list of datagrams, each has got five double-word (32-bit) values (DW0 to
                    # DW4). Little endian is storage is applied. Most significant byte (bits 31-24) of W0 tells the event type:
                    err = reinterpret(Int32, msg[1:4])[1]
                    overflow = reinterpret(UInt32, msg[5:8])[1]
                    usage = reinterpret(Float32, msg[9:12])[1]
                    encoding = reinterpret(UInt32, msg[13:16])[1]
                    ws = reinterpret(UInt32, msg[13:end])
                    frame = zero(eltype(ws))
                    for i = 1:length(ws)÷5
                        w0, w1, w2, w3, w4 = ws[5*(i-1)+1:5*i]
                        # println("$(ws[5*(i-1)+1:5*i])")
                        # msg[2] # timestamp in 10 ns units
                        msgid = (w0 >>> 24) & 0xFF
                        if msgid == 0x80 # dead time
                            # W0 = 0x80xxxxxx – dead-time measurement
                            # W0 – bits 0-7: filter index (1= filter 1, 2 = filter 2, …) bits 8-11: detector index (0=1st detector, 1=2nd detector, …)
                            # W1 – real-time timer value in 10 ns
                            # W2 – dead-time timer value in 10 ns
                            # W3 – event counter
                            x_det = 1 + ((w0 >>> 8) & 0x0F)
                            if dets[x_det] && frame in frames
                                @inbounds realtime += w1
                                @inbounds deadtime += w2
                            end
                        elseif msgid == 0x00 # measured energy
                            # W0 = 0x00xxxxxx – measured energy
                            # W0 – bits 0-7: filter index (1= filter 1, 2 = filter 2, …)  bits 8-11: detector index (0=1st detector, 1=2nd detector, ...)
                            # W1 – energy in eV (24.8 signed fixed point number)
                            # W2 – timestamp in 10ns
                            # W3 – pixel index
                            # W4 – frame index
                            x_det, frame = 1 + ((w0 >>> 8) & 0x0F), w4
                            if dets[x_det] && frame in frames
                                ch = channel(Float64(w1 >>> 8), scale)
                                r, c = Int(w3 ÷ rw - ry + 1)÷blocksize, Int(w3 % rw - rx + 1)÷blocksize
                                #elapsetime = ( prevelapse == typemax(UInt64) ? 0 : elapsetime + (w2 >= prevelapse ? UInt64(w2 - prevelapse) : zero(typeof(elapsetime))) # w2 + (prevelapse - typemax(UInt32)))
                                #prevelapse = w2
                                if checkbounds(Bool, hss.counts, ch, r, c)
                                    @inbounds hss.counts[ch, r, c] += 1
                                end
                            end
                        end
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
                    end
                elseif cmd == "ScData"
                    # See SharkSEM API ScScanLine documentation
                    # ScData(in unsigned int frameid, in int channel, in unsigned int index, in int bpp, in char[] data);
                    frameid = reinterpret(Int32, msg[1:4])[1]+1
                    maxframeid = max(frameid, maxframeid)
                    minframeid = min(frameid, minframeid)
                    ch = reinterpret(Int32, msg[5:8])[1]+1
                    maxchannel = max(ch, maxchannel)
                    index  = reinterpret(Int32, msg[9:12])[1]
                    bpp = reinterpret(Int32, msg[13:16])[1]
                    @assert bpp==8 || bpp==16
                    if !haskey(images,(frameid, ch))
                        images[(frameid, ch)] = (bpp == 8 ? zeros(UInt8, (rh - ry, rw - rx)) : zeros(UInt16, (rh - ry, rw - rx)))
                    end
                    img=images[(frameid, ch)]
                    imgdata = (bpp==8 ? reinterpret(UInt8, msg[21:end]) : reinterpret(UInt16, msg[21:end]))
                    img[index-ry+1:index-rx+length(imgdata)] .= imgdata
                else
                    if !(cmd in seen)
                        println("Skipping: $cmd")
                        push!(seen, cmd)
                    end
                end
                cmd, msg = readmsg(data)
            end
            if (maxchannel>0) && (maxframeid>0)
                for frameid in minframeid:maxframeid 
                    image = zeros(eltype(images[(minframeid,1)]),  (rh - ry, rw - rx, 1:maxchannel))
                    for ch in 1:maxchannel 
                        if haskey(images, (frameid, ch))
                            image[:,:,ch] .= images[(frameid,ch)]
                        end
                    end
                    hss[Symbol("Image$(maxframeid-minframeid+1)")] = gray.(permutedims(image, (2, 1, 3)))
                end
            end
            hss[:Dwell] = dwell * blocksize^2
            # Elapse time is broken!!!
            # hss[:Elapse] = 1.0e-8 * elapsetime # second
            hss[:RealTime] = (1.0e-8 * blocksize^2 / ((rh - ry) * (rw - rx))) * realtime  # seconds
            hss.livetime .= (1.0e-8 * blocksize^2 / ((rh - ry) * (rw - rx))) * (realtime - deadtime)  # seconds
            hss[:DeadFraction] = deadtime / realtime
            return hss
        end
    end # Delete the temporary directory
end
