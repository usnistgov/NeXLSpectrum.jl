using CSV

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
table.

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

function readrpl(io::IO)::RPLHeader
    w, h, d, o = -1, -1, -1, -1
    dl, dtv, bo, rb = -1, "", :unknown, :unknown
    for row in CSV.Rows(read(io))
        key, val = uppercase(strip(row[1])), strip(row[2])
        if key == "WIDTH"
            w = parse(Int, val)
        elseif key == "HEIGHT"
            h = parse(Int, val)
        elseif key == "DEPTH"
            d = parse(Int, val)
        elseif key == "OFFSET"
            o = parse(Int, val)
        elseif key == "DATA-LENGTH"
            dl = parse(Int, val)
        elseif key == "DATA-TYPE"
            dtv = uppercase(val)
        elseif key == "BYTE-ORDER"
            bo = uppercase(val) == "BIG-ENDIAN" ? :bigendian : :littleendian
        elseif key == "RECORD-BY"
            rb = uppercase(val) == "IMAGE" ? :image : :vector # :dontcare => :vector
        end
    end
    dt = UInt8
    if dtv == "SIGNED"
        if dl == 1
            dt = Int8
        elseif dl == 2
            dt = Int16
        elseif dl == 4
            dt = Int32
        elseif dl == 8
            dt = Int64
        else
            error("Unexpected data-length $dl")
        end
    elseif dtv == "UNSIGNED"
        if dl == 1
            dt = UInt8
        elseif dl == 2
            dt = UInt16
        elseif dl == 4
            dt = UInt32
        elseif dl == 8
            dt = UInt64
        else
            error("Unexpected data-length $dl")
        end
    elseif dtv == "FLOAT"
        if dl == 2
            dt = Float16
        elseif dl == 4
            dt = Float32
        elseif dl == 8
            dt = Float64
        else
            error("Unexpected data-length $dl")
        end
    else
        error("Unexepected data-type $dt.")
    end
    return RPLHeader(w, h, d, o, dt, bo, rb)
end

"""
    readraw(ios::IOStream, rpl::RPLHeader)::Array{<:Real}

Construct an Array from the data in binary from the specified stream.  Always returns the data in [d,w,h] order
regardless of how the data is organized on disk.
"""
function readraw(ios::IOStream, rpl::RPLHeader)::Array{<:Real}
    if rpl.recordedBy == :vector
        size = ( rpl.depth, rpl.height, rpl.width )
    else
        size = ( rpl.height, rpl.width, rpl.depth )
    end
    seek(ios, rpl.offset)
    data = Array{rpl.dataType}(undef, size...)
    # Data in a RPL file is organized ch -> row -> col
    for r in Base.OneTo(rpl.height), c in Base.OneTo(rpl.width) 
        read!(ios, @view data[:, r, c])
    end
    nativeendian = (ENDIAN_BOM == 0x04030201 ? :littleendian : :bigendian)
    if rpl.byteOrder != nativeendian
        if nativeendian == :littleendian
            foreach(i -> data[i] = ntoh(data[i]), eachindex(data))
        else
            foreach(i -> data[i] = hton(data[i]), eachindex(data))
        end
    end
    res = data
    if rpl.recordedBy == :image
        # Reorder as :vector
        res = Array{rpl.dataType}(undef, rpl.depth, rpl.width, rpl.height)
        for d in 1:rpl.depth, h in 1:rpl.height, w in 1:rpl.width 
            res[d, h, w] = data[h, w, d]
        end
    end
    return res
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
