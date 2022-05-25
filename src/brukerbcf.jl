# Read the Bruker AidAim SingleFileSystem based BCF hyperspectrum files.
# Format was partially reverse engineered (partly - as there are still some
# not well defined parts of the format) and below parsing implementation was
# written by Petras Jokubauskas,
# who dedicates this implementation to the public domain (Unlicense).
using Dates

const SFS_MAGIC_STRING = Vector{UInt8}("AAMVHFSS")


"""
    readbrukerbcf(filename::AbstractString,
                  getelectronimage::Bool,
                  downsample::UInt8)
                  ::[Array{UInt16, (w::UInt16, h::UInt16)}, HyperSpectra]
    by-default: getelectronimage is true
                downsample == 1 (no downsampling)

Read a Bruker BCF EDS hyperspectrum file into a `HyperSpectrum` and optionaly
return the Scanning Electron Image(s) (BSE, SEI, In-Lens, Argus...
depends from which non-Xray detector(s) was/were selected for acquisition).
"""
function readbrukerbcf(filename::String, getelectronimage::Bool=true, downsample::UInt8=1)::Spectrum
    return open(filename, read = true) do io
        readbrukerspx(io)
    end
end

function msfiletimetodatetime(filetime::UInt64)::DateTime
    unix_ts = filetime / 10000000 - 11644473600
    return unix2datetime(unix_ts)
end

"""SFS Header is a container of read and precalculated values which later
are needed for reading any embedded file. The nomenclature here uses chunk
probably page would be more correct. chunk indexes starts with 0, at 0th chunk
the SFS header is situated, which describes chunk size, index of chunk with
file tree table, total number of chunks.
This or part of information is required for any reading.
chunsize are forsed to UInt64, so the later resolution into full pointer would not overflow.
"""
struct SfsHeader
    version::String
    chunksize::UInt64  # 
    chunkdatasize::UInt64 # chunksize minus header 
    sfs_tree_first_chunk:: UInt64 # # chunk number 
    n_filetree_items::UInt32
    n_sfs_chunks::UInt32
end

"""read the Aid Aim Single File System file and return hierarchical dictionary with files and chunk address tables neccesarry for reading the embedded files"""
function readsfs(filename::String)
    open(filename, "r") do io
        file_signature = read(io, 8)
        if file_signature ≠ SFS_MAGIC_STRING
            throw(ArgumentError("Not a Bruker bcf file, file $filename does not begin with 8byte signature 0xAAMVHFSS (Aid Aim Single FileSystem)"))
        end
        # FEFFFFFF FFFFFF7F - meaning of that is unknown
        seek(io, 0x124)  # up to 0x118 is zeroes; 0x118 starts with another signature 
        version = read(io, Float32)
        chunksize = read(io, UInt32)  # depends form versions of sfs and Esprit
        usable_chunksize = chunksize - 0x20
        seek(io, 0x140)
        # location of sfs hierarchical tree first chunk, number of
        # items in the tree, and total number of chunks in the sfs file
        tree_start_chunk = read(io, UInt32) # the index of the chunk
        n_tree_items = read(io, UInt32)
        n_of_sfs_chunks = read(io, UInt32)
        sfsheader = SfsHeader(string(version), chunksize,
                              usable_chunksize, tree_start_chunk,
                              n_tree_items, n_of_sfs_chunks)
        n_tree_chunks = ceil(UInt32,(n_tree_items * 0x200) / (chunksize - 0x20))
        if n_tree_chunks == 1 # file tree do not exceed one chunk in bcf
            seek(io, chunksize * tree_start_chunk + 0x138)
            tree_buffer = IOBuffer(read(io, 0x200 * n_tree_items))
        # in simple Hypermap situation file tree fits into single block
        # however, in case we are intereseted in hypermaps gathered together
        # with EBSD, such file tree would extend few blocks as EBSD produces
        # dozens of files
        else
            tree_buffer = IOBuffer(read=true, write=true, append=true)
            tree_pointer = tree_start_chunk
            tree_items_in_chunk = floor(UInt32, usable_chunksize / 0x200)
            for _=1:tree_items_in_chunk
                # jump to chunk header:
                seek(io, chunksize * tree_pointer + 0x118)
                #set the next chunk pointer:
                tree_pointer = read(io, UInt32)
                skip(io, 28) # probably some crt or checksums, not well defined
                write(tree_buffer, read(io, tree_items_in_chunk * 0x200))
            end
            seekstart(tree_buffer)
        end
        files = Vector{SfsItem}()
        for _=1:n_tree_items
            push!(files, SfsItem(tree_buffer, io, sfsheader))
        end
        return sfsheader, files
    end
end

struct SfsItem
    filesize::UInt64  # bytes the file is taking on disk
    createtime::DateTime
    modificationtime::DateTime
    accesstime::DateTime
    fileflags::UInt32
    parent::Int32 # index pointer (zero based! -1 points to no parent)
    isdir::Bool
    name::String
    pointertable::Vector{UInt64}  # internaly sfs pointer tables are 32bit
    # final location of chunk is calulated by multiplying with chunks size
    # and adding an fixed offset of 0x118; here the pointertable is 
    # with ready final pointers in 64bits (BCF files can far exceed 3.5GB
    # limitation of 32bit addressing resolution)
    function SfsItem(iotree::IO, iomain::IO, header::SfsHeader)
        startpos = position(iotree)
        tablepointer = UInt64(read(iotree, UInt32))
        filesize = read(iotree, UInt64)
        ctime = msfiletimetodatetime(read(iotree, UInt64))
        mtime = msfiletimetodatetime(read(iotree, UInt64))
        atime = msfiletimetodatetime(read(iotree, UInt64))
        flags = read(iotree, UInt32)  # read, write, execute permissions
        parent = read(iotree, Int32)
        skip(iotree, 176)  # encryption stuff? normally for bcfs all 0
        isdir = read(iotree, Bool)
        skip(iotree, 3)  # bool is written over 4 bytes (4-1=3)
        name = readuntil(iotree, "\0") # Null byte terminated C str (=<256)
        # meaning of last 32 bytes are not known
        n_chunks = ceil(UInt32, filesize / header.chunkdatasize)
        file_ptr_table = UInt64[]
        if !isdir
            # table size in number of chunks:
            n_tab_chunks = ceil(UInt32, n_chunks * 4 / header.chunkdatasize)
            tab_iobuffer = IOBuffer(read=true, write=true, append=true)
            if n_tab_chunks > 1
                next_pointer = tablepointer
                for _=1:n_tab_chunks
                    seek(iomain, next_pointer * header.chunksize + 0x118)
                    next_pointer = UInt64(read(iomain, UInt32))
                    skip(iomain, 28) # skip rest of the chunk header
                    write(tab_iobuffer, read(iomain, header.chunkdatasize))
                end
            else
                seek(iomain, tablepointer * header.chunksize + 0x138)
                # 0x118 general offset + 0x20 (32 bytes) from skipping
                # whole header part
                write(tab_iobuffer,(read(iomain, header.chunkdatasize)))

            end
            seekstart(tab_iobuffer)
            for _=1:n_chunks
                push!(file_ptr_table,
                      UInt64(read(tab_iobuffer, UInt32)) * header.chunksize + 0x118)
            end     
        end
        seek(iotree, startpos + 0x200) # prepare IO for other item read
        new(filesize, ctime, mtime, atime,
            flags, parent, isdir, name,
            file_ptr_table)
    end
end
