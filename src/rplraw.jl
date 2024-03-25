using Unitful

struct RPLHeader
    width::Int
    height::Int
    depth::Int
    offset::Int
    dataType::Type{<:Real}
    byteOrder::Symbol # :bigendian, :littleendian
    recordedBy::Symbol # :vector, :image
    function RPLHeader(
        w::Int,
        h::Int,
        d::Int,
        o::Int,
        dt::Type{<:Real},
        bo::Symbol,
        rb::Symbol,
    )
        @assert w > 0
        @assert h > 0
        @assert d > 0
        @assert o >= 0
        @assert dt in [
            Int8,
            Int16,
            Int32,
            Int64,
            UInt8,
            UInt16,
            UInt32,
            UInt64,
            Float16,
            Float32,
            Float64,
        ]
        @assert bo in [:littleendian, :bigendian]
        @assert rb in [:vector, :image]
        new(w, h, d, o, dt, bo, rb)
    end
end

readrplraw(rplio::IO, rawio::IO) = readraw(rawio, readrpl(rplio))

"""
    readrplraw(rplfilename::AbstractString, rawfilename::AbstractString)
    readrplraw(rplio::IO, rawio::IO)

Read a RPL/RAW file pair from IO or filename into an Array obect.  The reader supports :bigendian, :littleendian
ordering and :vector or :image alignment.  Since the Array can be very large it is beneficial to release and collected
the associated memory when you are done with the data by assigning `nothing` to the variable (or allowing it to go
out of scope) and then calling `GC.gc()` to force a garbage collection.


    * Ordering: The individual data items may be in `:littleendian` or `:bigendian`.
    ** `:littleendian` => Intel/AMD and others
    ** `:bigendian` => ARM/PowerPC/Motorola
    * Alignment:  The data in the file can be organized as `:vector` or `:image`.  However, the data will be
reorganized into 'vector' format when returned as a Array.
    **`:vector` => Spectrum/data vector as contiguous blocks by pixel
    ** `:image` => Each channel of data organized in image planes
    * Data types: signed/unsigned 8/16/32-bit integers or 16-bit/32-bit/64-bit floats

    readrplraw(filenamebase::AbstractString)::Array{<:Real}

Read the files filenamebase*".rpl" and filenamebase*".raw" into an Array.  Maintains the data type
of the values in the RAW file.


#### Standard LISPIX Parameters in .rpl File

.rpl files consist of 9 lines.  Each line consists of a 'key'<tab>'value' where there is one and only one tab and
possibly other space between the parameter name and parameter value. Parameter names are case-insensitive.
The first line in the files is "key<tab>value".  Subsequent lines contain the keys and values described in this
table.  Other programs (namely HyperSpy) have stuffed additional data into the RPL header using a similar syntax.
HyperSpy also uses ';' to start comment lines.

| **key**      |  **value**    | **description**                          |
| :------------|---------------|:-----------------------------------------|
| width        | 849           |  pixels per row       integer            |
| height       | 846           |  rows                 integer            |
| depth        | 4096          |  images or spec pts   integer            |
| offset       | 0             |  bytes to skip        integer            |
| data-length  | 1             |  bytes per pixel      1, 2, 4, or 8      |
| data-type    | unsigned      |  signed, unsigned, or float              |
| byte-order   | dont-care     |  big-endian, little-endian, or dont-care |
| record-by    | vector        |  image, vector, or dont-care             |

This .rpl file indicates the image is 849 pixels wide and 846 pixels high, with 4096 levels in the depth dimension.

Example HyperSpy RPL file:

```
key	value
beam-energy	12.0
byte-order	dont-care
data-length	1
data-type	unsigned
date	2024-01-11
depth	4096
depth-name	Energy
depth-origin	-0.47665
depth-scale	0.005
depth-units	keV
elevation-angle	35.0
ev-per-chan	0
height	1024
height-name	height
height-origin	0.0
height-scale	0.9734798568
height-units	µm
offset	0
record-by	vector
signal	EDS_SEM
tilt-stage	-0.001
time	12:21:50
title	EDX
width	1024
width-name	width
width-origin	0.0
width-scale	0.9734798568
width-units	µm
```
"""
function readrplraw(rplfilename::AbstractString, rawfilename::AbstractString)::Array{<:Real}
    open(rplfilename, read = true) do rplio
        open(rawfilename, read = true) do rawio
            return readrplraw(rplio, rawio)
        end
    end
end

readrplraw(filenamebase::AbstractString)::Array{<:Real} =
    readrplraw(String(filenamebase) * ".rpl", String(filenamebase) * ".raw")


"""
Reads either Bright- or Oxford-style RPL headers.
"""
function readrpl(io::IO)::RPLHeader
    w, h, d, o = -1, -1, -1, -1
    dl, dtv, bo, rb, dt = -1, "", :unknown, :unknown, Nothing
    items = Dict{String,String}()
    for line in split(read(io, String),"\n")
        kv = if startswith(line,"(MLX:") # Oxford-style
            split(line[7:end-2]," ")
        else 
            split(line,"\t") # Bright-style
        end
        if length(kv)!=2
            continue
        end
        key, val = uppercase(kv[1]), strip(kv[2])
        items[key] = val
    end
    if haskey(items, "WIDTH")
        w = parse(Int, items["WIDTH"])
    end
    if haskey(items, "HEIGHT")
        h = parse(Int, items["HEIGHT"])
    end
    if haskey(items, "DEPTH")
        d = parse(Int, items["DEPTH"])
    end
    if haskey(items, "OFFSET")
       o = parse(Int, items["OFFSET"])
    end
    if haskey(items, "DATA-LENGTH")
       dl = parse(Int, items["DATA-LENGTH"])
    end
    if haskey(items, "DATA-TYPE")
       dtv = uppercase(items["DATA-TYPE"])
       dt = if dtv in ( "SIGNED", ":SIGNED" )
           if dl == 1
               Int8
           elseif dl == 2
               Int16
           elseif dl == 4
               Int32
           elseif dl == 8
               Int64
           else
               Nothing
           end
       elseif dtv in ( "UNSIGNED", ":UNSIGNED" )
           if dl == 1
               UInt8
           elseif dl == 2
               UInt16
           elseif dl == 4
               UInt32
           elseif dl == 8
               UInt64
           else
               Nothing
           end
       elseif dtv in ( "FLOAT", ":FLOAT" )
           if dl == 2
               Float16
           elseif dl == 4
               Float32
           elseif dl == 8
               Float64
           else
               Nothing
           end
       else
           Nothing
       end
    end
    if dt == Nothing
        error("Unexpected data-length $dl")
    end
    if haskey(items, "BYTE-ORDER")
        bo = uppercase(items["BYTE-ORDER"]) in ( "BIG-ENDIAN", ":BIG-ENDIAN") ? :bigendian : :littleendian
    end
    if haskey(items, "RECORD-BY")
        rb = uppercase(items["RECORD-BY"]) in ( "IMAGE", ":IMAGE") ? :image : :vector # :dontcare => :vector
    end
    return RPLHeader(w, h, d, o, dt, bo, rb)
end

"""
    readraw(ios::IOStream, rpl::RPLHeader)::Array{<:Real}

Construct an Array from the data in binary from the specified stream.  Always returns the data in [d,w,h] order
regardless of how the data is organized on disk.
"""
function readraw(ios::IOStream, rpl::RPLHeader)::Array{<:Real}
    seek(ios, rpl.offset)
    data = Array{rpl.dataType}(undef, rpl.depth*rpl.height*rpl.width)
    # Data in a RPL file is organized ch -> row -> col or row -> col -> ch
    read!(ios, data)
    nativeendian = (ENDIAN_BOM == 0x04030201 ? :littleendian : :bigendian)
    if rpl.byteOrder != nativeendian
        if nativeendian == :littleendian
            map!(ntoh, data, data)
        else
            map!(hton, data, data)
        end
    end
    if rpl.recordedBy == :image
        # Reorder as :vector
        vector = permutedims(reshape(data, (rpl.width, rpl.height, rpl.depth)), (3, 2, 1))
    else 
        vector = reshape(data, (rpl.depth, rpl.height, rpl.width))
    end
    return vector
end

function writerplraw(rplbasefile::String, arr::AbstractArray{<:Real})
    open(rplbasefile * ".rpl", write = true) do rplio
        println(rplio, "key\tvalue")
        println(rplio, "width\t$(size(arr,2))")
        println(rplio, "height\t$(size(arr,3))")
        println(rplio, "depth\t$(size(arr,1))")
        println(rplio, "offset\t0")
        println(rplio, "data-length\t$(sizeof(eltype(arr)))")
        if eltype(arr) in (Int8, Int16, Int32, Int64)
            println(rplio, "data-type\tsigned")
        elseif eltype(arr) in (UInt8, UInt16, UInt32, UInt64)
            println(rplio, "data-type\tunsigned")
        elseif eltype(arr) in (Float16, Float32, Float64)
            println(rplio, "data-type\tfloat")
        else
            error("Unsupported type $(eltype(arr)) in write RPL/RAW.")
        end
        println(rplio, "byte-order\tlittle-endian")
        println(rplio, "record-by\tvector")
    end
    open(rplbasefile * ".raw", write = true) do rawio
        write(rawio, arr)
    end
end

function writerplraw(rplbasefile::String, hs::HyperSpectrum{T,2,3}) where { T <: Real }
    # Understand length and time units
    to_μm(sc::Quantity) =   dimension(sc)==dimension(u"m") ? ( ustrip(uconvert(u"μm",sc)), "μm" ) : ( dimension(sc)==dimension(u"s") ? ( ustrip(uconvert(u"s",sc)), "s" ) :  ( ustrip(sc), "unknown" ) )
    to_μm(sc::Real) = ( sc, "unknown" ) 
    open(rplbasefile * ".rpl", write = true) do rplio
        println(rplio, "key\tvalue")
        e0 = get(hs, :BeamEnergy, nothing)
        if !isnothing(e0)
            println(rplio, "beam-energy\t$(e0/1.0e3)")
        end
        println(rplio, "width\t$(size(hs,1))")
        println(rplio, "width-name\twidth")
        ( orig, unt ) = to_μm(axisvalue(hs, 2, 1))        
        ( sc, _  ) = to_μm((axisvalue(hs, 2, 2)-axisvalue(hs, 2, 1)))
        println(rplio, "width-origin\t$orig")
        println(rplio, "width-scale\t$sc")
        println(rplio, "width-units\t$unt")
        println(rplio, "height\t$(size(hs,2))")
        println(rplio, "height-name\theight")
        ( orig, unt ) = to_μm(axisvalue(hs, 2, 1))        
        ( sc, _  ) = to_μm((axisvalue(hs, 2, 2)-axisvalue(hs, 2, 1)))
        println(rplio, "height-origin\t$orig")
        println(rplio, "height-scale\t$sc")
        println(rplio, "height-units\t$unt")
        println(rplio, "depth\t$(depth(hs))")
        println(rplio, "depth-name\tEnergy")
        println(rplio, "depth-origin\t$(energy(1,hs)/1000.0)")
        println(rplio, "depth-scale\t$(channelwidth(1,hs)/1000.0)")
        println(rplio, "depth-units\tkeV")
        println(rplio, "offset\t0")
        println(rplio, "data-length\t$(sizeof(T))")
        if T in (Int8, Int16, Int32, Int64)
            println(rplio, "data-type\tsigned")
        elseif T in (UInt8, UInt16, UInt32, UInt64)
            println(rplio, "data-type\tunsigned")
        elseif T in (Float16, Float32, Float64)
            println(rplio, "data-type\tfloat")
        else
            error("Unsupported type $T in write RPL/RAW.")
        end
        println(rplio, "byte-order\tlittle-endian")
        println(rplio, "record-by\tvector")
        θ = get(hs, :TakeOffAngle, nothing)
        if !isnothing(θ)
            println(rplio, "elevation-angle\t$(rad2deg(θ))")
        end
        acq = get(hs,:AcquisitionTime,nothing)
        if !isnothing(acq)
            dt = "$(Dates.year(acq))-$(Dates.month(acq))-$(Dates.day(acq))"
            tm = "$(Dates.hour(acq)):$(Dates.minute(acq)):$(Dates.second(acq))"
            println(rplio, "date\t$dt")
            println(rplio, "time\t$tm")
        end
        println(rplio, "signal\tEDS_SEM")
        # println(rplio, "tilt-stage	-0.001
        title = get(hs, :Title, "EDX")
        println(rplio, "title\t$title")
        println(rplio, "offset\t0")
        println(rplio, "ev-per-chan\t0")
    end
    open(rplbasefile * ".raw", write = true) do rawio
        write(rawio, hs.counts)
    end
end