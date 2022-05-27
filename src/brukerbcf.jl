# Read the Bruker AidAim SingleFileSystem based BCF hyperspectrum files.
# Format was partially reverse engineered (partly - as there are still some
# not well defined parts of the format) and below parsing implementation was
# written by Petras Jokubauskas,
# who dedicates this implementation to the public domain (Unlicense).
using Dates 

####################################SFS part####################################
# Bruker element hypermaps (.bcf) (bcf is  lso used for EBSD) and
# particle analysis files (.pan) are built on top of AidAim Software virtual 
# SingleFileSystem (further called sfs). It allows overcome 32 bit limitiation 
# (Esprit1.x 32bit version) for file size, and most importantly that technology
# allows memory efficient random access to embedded files,
# despite using of compression (i.e. accessing a single pixel spectrum without
# loading whole 3D array into memory).
#
# The implementation of SFS is only partial and simplified - it contains just
# a subset of functionality mapping to SFS features used by Bruker files.
# As .bcf (and .pan) use shallow hierarchy, and there is known name collisions
# - no hierarchical reconstruction of sfs file tree is implemented
# (File pointer tables are returned in 1D vector).
#
# SFS divides the file into `chunks` and that is primary blocks of SFS
# (also could be called "Pages" in file system context). Initial (0th) block
# contains necessary information about chunk size, total number of chunks,
# version, index of chunk with start offile tree. Every chunk starts with
# 32 byte chunk header which meaning is not fully reverse engineered (probably
# contains some checksums).
#
# File tree can span multiple chunks - the items (SfsItem) contain parent index,
# is dir flag, chunk index where file pointer table starts, file/folder name.
# Pointers are chunk indices, real pointer for reading content needs to be
# made by multiplying that with chunksize, and offset with 0x118 (universal
# offset of all chunks).
#
# In case of compression, Esprit uses only zlib compression (although bz2 and 
# other compression is possible in SFS) and thus only two cases: not compressed 
# and zlib compressed data are handled. Compression in SFS works as another
# layer which divides files into the blocks (that allows to extract only 
# needed chunk of data at the middle of the data stream).
# That however complicates the reading of file contents, as blocks of
# compression layer can be crossing chunk boundaries of SFS file.
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

const SFS_MAGIC_STRING = Vector{UInt8}("AAMVHFSS")

"""Convert Win32 FILETIME to Julia DateTime by converting the FILETIME
into UNIX timestamp and then converting that into datetime"""
function msfiletimetodatetime(filetime::UInt64)::DateTime
    unix_ts = filetime / 10000000 - 11644473600
    return unix2datetime(unix_ts)
end

"""SFS Header is a container of read and precalculated values which later
are needed for reading any embedded file.
The chunk indexes starts from 0. Values are filled into this struct from that
part, as header is small enought to fit into most tiny chunk size imaginable.
Recognised fields are:
- version (good to know for debuggin in case of problems in future),
- chunksize,
- sfs_tree_first_chunk - chunk index where file tree table begins,
- n_filetree_items - number of file tree items,
- n_sfs_chunks -total number of chunks.

Additionally struct contains calculated and inhertied values
chunkdatasize - size without chunk header.
file_path - for reopening file

This or part of information is required for any reading.
chunk(data)size(s) are forced to UInt64, so the later resolution into
full address pointer would not overflow 32bits.
"""
struct SfsHeader
    version::Float32
    chunksize::UInt64   
    chunkdatasize::UInt64 # chunksize minus header 
    sfs_tree_first_chunk:: UInt64 # chunk index 
    n_filetree_items::UInt32
    n_sfs_chunks::UInt32
    filename::String
end

"""parse the Aid Aim Single File System file and return representation of
single file system as SFVFS (single file virtual file system) type,
that contains necessary information for extracting the content
of the first layer of SFS
"""
function getsfvfs(filename::String)::SFVFS
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
        sfsheader = SfsHeader(version, chunksize,
                              usable_chunksize, tree_start_chunk,
                              n_tree_items, n_of_sfs_chunks,
                              filename)
        n_tree_chunks = cld(n_tree_items * 0x200, usable_chunksize)
        if n_tree_chunks == 1 # file tree do not exceed one chunk in bcf
            seek(io, chunksize * tree_start_chunk + 0x138)
            tree_buffer = IOBuffer(read(io, 0x200 * n_tree_items))
        # in simple bcf Hypermap situation whole file tree fits into 1 block.
        # However, in case we are intereseted in hypermaps gathered together
        # with EBSD, such file tree would extend few blocks as EBSD produces
        # dozens of files. Same SFS technology is also used for .pan
        # (particle analysis) files, and there file table could easily span
        # dozens of chunks
        else
            tree_buffer = IOBuffer(read=true, write=true, append=true)
            tree_pointer = tree_start_chunk
            tree_items_in_chunk = fld(usable_chunksize, 0x200)
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
        return SFVFS(sfsheader, files)
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
    pointertable::Vector{UInt64}  # internaly sfs pointer tables are 32bit,
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
        n_chunks = cld(filesize, header.chunkdatasize)
        file_ptr_table = UInt64[]
        if !isdir
            # table size in number of chunks:
            n_tab_chunks = cld(n_chunks * 4, header.chunkdatasize)
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

"""representation of single file virtual file system"""
struct SFVFS
    header::SfsHeader
    contents::Vector{SfsItem}
end

function getsfsitem(list::Vector{SfsItem}, name::String)::SfsItem
    for i in list
        if name == i.name
            return i
        end
    end
end

"""
returns the selected by name file from SFVFS as IOBuffer
"""
function sfsitem2iobuffer(sfs::SFVFS, name::String)::IO
   item = getsfsitem(sfs.contents, name)
   if !isnothing(item)
       io = IOBuffer()
       remaining_bytes = item.filesize
       open(sfs.header.filename, "r") do fh
           for chunk_ptr in item.pointertable
               if remaining_bytes > sfs.header.chunkdatasize
                   n_to_read = sfs.header.chunkdatasize
               else
                   n_to_read = remaining_bytes
               end
               seek(fh, chunk_ptr + 0x20)
               write(io, read(fh, n_to_read))
               remaining_bytes -= n_to_read
           end
       end
       seek(io, 0)
       return io
   end
end


mutable struct SFSItemIO <: IO
    item::SfsItem
    opened::Bool
    pos::UInt64
    size::UInt64
    innerio::IO
    function SFSItemIO(sfs::SFVFS, name::String)
        item = getsfsitem(sfs.contents, name)
        if !isnothing(item)
            io = open(sfs.header.filename, "r")
            new(item, true, UInt64(0), io)
        end
    end
end

######################################BCF part#################################
# Bcf uses the SFS as file container.
# Element Hypermaps are organised inside the SFS like this:
# / (root_folder)
# └─EDSDatabase
#   ├─HeaderData
#   ├─SpectrumData0
#   └─SpectrumPositions0     (Starting with Esprit version 2)
#
# In case of Elemental Hypermaps embedded together with EBSD in a single bcf
# the structure will look like this:
# / (root_folder)
# ├─EDSDatabase
# │ ├─HeaderData
# │ ├─SpectrumData0
# │ └─SpectrumPositions0 
# └─EBSDData
#   ├─SEMImage
#   ├─SEMStageData
#   ├─Calibration
#   ├─Subsets
#   ├─GrainDetection
#   ├─MapMask
#   ...(30 and more files)
# 
# HeaderData - xml file with embedded sum spectra of whole Hyperspectra
# and image nodes (containing SEM TEM electron (or ect) 16bit image(s))
# SpectrumData0 - binary format containing the HyperMap saved as
# Packed Delphi Pascal Array with internal offset of 0x1A0. "Packed"
# should be very much emphasized as it makes it impossible to simply memmap,
# but requires meticulous unpacking (including lots of bit shiftting).
# SpectrumPositions0 - due to very dynamic nature of packed Delphi arrays,
# getting to particular (x,y) pixel would require to parse spectral pixel from
# (0,0). SpectrumPositions0 alows to access pixel spectrum without parsing
# the whole HyperMap. It contains array of pixel pointers
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV 

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


