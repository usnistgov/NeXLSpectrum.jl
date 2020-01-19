using CSV

struct RPLHeader
    width::Int
    height::Int
    depth::Int
    offset::Int
    dataType::Type{<:Real}
    byteOrder::Symbol # :bigendian, :littleendian
    recordedBy::Symbol # :vector, :image
    function RPLHeader(w::Int,h::Int,d::Int,o::Int,dt::Type{<:Real},bo::Symbol,rb::Symbol)
        @assert w>0
        @assert h>0
        @assert d>0
        @assert o>=0
        @assert dt in [ Int8, Int16, Int32, UInt8, UInt16, UInt32, Float16, Float32, Float64 ]
        @assert bo in [ :littleendian, :bigendian ]
        @assert rb in [ :vector, :image ]
        new(w,h,d,o,dt,bo,rb)
    end
end



readrplraw(rplio::IO, rawio::IO, scale::EnergyScale, props::Dict{Symbol,Any}=Dict{Symbol,Any}()) =
    readraw(rawio, readrpl(rplio), scale, props)

"""
    readrplraw(rpl, raw, scale::EnergyScale, props::Dict{Symbol,Any})

Read a RPL/RAW file pair from IO or File (by name).
"""
function readrplraw(rplfilename::AbstractString, rawfilename::AbstractString, scale::EnergyScale, props::Dict{Symbol,Any}=Dict{Symbol,Any}())
    open(rplfilename, read=true) do rplio
        open(rawfilename, read=true) do rawio
            return readrplraw(rplio, rawio, scale, props)
        end
    end
end

readrplraw(filenamebase::AbstractString, scale::EnergyScale, props::Dict{Symbol,Any}=Dict{Symbol,Any}()) =
    readrplraw(String(filenamebase)*".rpl", String(filenamebase)*".raw", scale, props)

function readrpl(io::IO)
    w, h, d, o = -1, -1, -1, -1
    dl, dtv, bo, rb = -1, "", :unknown, :unknown
    for row in CSV.File(io)
        key, val = uppercase(strip(row[1])), strip(row[2])
        if key=="WIDTH"
            w=parse(Int,val)
        elseif key=="HEIGHT"
            h=parse(Int,val)
        elseif key=="DEPTH"
            d=parse(Int,val)
        elseif key=="OFFSET"
            o=parse(Int,val)
        elseif key=="DATA-LENGTH"
            dl=parse(Int,val)
        elseif key=="DATA-TYPE"
            dtv=val
        elseif key=="BYTE-ORDER"
            bo = uppercase(val)=="BIG-ENDIAN" ? :bigendian : :littleendian
        elseif key=="RECORD-BY"
            rb = uppercase(val)=="IMAGE" ? :image : :vector # :dontcare => :vector
        end
    end
    dt=UInt8
    if dtv=="SIGNED"
        if dl == 1
            dt=Int8
        elseif dl == 2
            dt=Int16
        elseif dl == 4
            dt=Int32
        elseif dl==8
            dt=Int64
        else
            error("Unexpected data-length $dl in $rplfilename")
        end
    elseif dtv=="UNSIGNED"
        if dl == 1
            dt=UInt8
        elseif dl == 2
            dt=UInt16
        elseif dl == 4
            dt=UInt32
        elseif dl==8
            dt=UInt64
        else
            error("Unexpected data-length $dl in $rplfilename")
        end
    elseif dtv=="FLOAT"
        if dl==2
            dt=Float16
        elseif dl == 4
            dt=Float32
        elseif dl==8
            dt=Float64
        else
            error("Unexpected data-length $dl in $rplfilename")
        end
    else
        error("Unexepected data-type $dt in $rplfilename.")
    end
    return RPLHeader(w,h,d,o,dt,bo,rb)
end

"""
    readraw(ios::IOStream, T::Type{<:Real}, size::NTuple{Int, N}, energy::EnergyScale, props::Dict{Symbol,Any}) =

Construct a full Signal from the counts data in binary from the specified stream.  Always
returns the data in [d,w,h] order.
"""
function readraw(ios::IOStream, rpl::RPLHeader, energy::EnergyScale, props::Dict{Symbol,Any})
    if rpl.recordedBy==:vector
        size = ( rpl.depth, rpl.width, rpl.height )
    else
        size = ( rpl.width, rpl.height, rpl.depth )
    end
    seek(ios,rpl.offset)
    data = read!(ios, Array{rpl.dataType}(undef,size...))
    nativeendian = (ENDIAN_BOM==0x04030201 ? :littleendian : :bigendian)
    if rpl.byteOrder != nativeendian
        if nativeendian==:littleendian
            foreach(i->data[i]=ntoh(data[i]),eachindex(data))
        else
            foreach(i->data[i]=hton(data[i]),eachindex(data))
        end
    end
    res = data
    if rpl.recordedBy==:image
        # Reorder as :vector
        res = Array{rpl.dataType}(undef, rpl.depth, rpl.width, rpl.height)
        for d in 1:rpl.depth, w in 1:rpl.width, h in 1:rpl.height
            res[d,w,h] = data[w,h,d]
        end
    end
    return Signal(energy, props, res)
end

function writerplraw(rplbasefile::String, sig::Signal)
    open(rplbasefile*".rpl",write=true) do rplio
        println(rplio,"key\tvalue")
        println(rplio,"width\t$(size(sig,2))")
        println(rplio,"height\t$(size(sig,3))")
        println(rplio,"depth\t$(size(sig,1))")
        println(rplio,"offset\t0")
        println(rplio,"data-length\t$(sizeof(sig.counts[1]))")
        if typeof(sig.counts[1]) in (Int8, Int16, Int32)
            println(rplio,"data-type\tsigned")
        elseif typeof(sig.counts[1]) in (UInt8, UInt16, UInt32)
            println(rplio,"data-type\tunsigned")
        elseif typeof(sig.counts[1]) in (Float16, Float32, Float64)
            println(rplio,"data-type\tfloat")
        else
            error("Unexpected type $(typeof(sig.counts[1])) in write RPL/RAW.")
        end
        println(rplio,"byte-order\tlittle-endian")
        println(rplio,"record-by\tvector")
    end
    open(rplbasefile*".raw",write=true) do rawio
        write(rawio, sig.counts)
    end
end
