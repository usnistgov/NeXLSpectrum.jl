# Read the Bruker AidAim SingleFileSystem based BCF hyperspectrum files.
# Format was partially reverse engineered (partly - as there are still some
# not well defined parts of the format) and below parsing implementation was
# written by Petras Jokubauskas,
# who dedicates this implementation to the public domain (Unlicense).

const SFS_HEADER = Vector{UInt8}("AAMVHFSS")


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


"""read the Aid Aim Single File System file and return hierarchical dictionary with files and chunk address tables neccesarry for reading the embedded files"""
function readsfs(filename::String)
    open(filename, "r") do io
        file_signature = read(io, 8)
        if file_signature ≠ SFS_HEADER
            throw(ArgumentError("Not a Bruker bcf file, file $filename does not begin with 8byte signature 0xAAMVHFSS (Aid Aim Single FileSystem)"))
        end
        seek(io, 0x124)
        version = read(io, Float32)
        chunksize = read(io, UInt32)
        usable_chunksize = chunksize - 0x20
        seek(io, 0x140)
        # location of sfs hierarchical tree first chunk, number of
        # items in the tree, and total number of chunks in the sfs file
        tree_start_chunk = read(io, UInt32) # the index of the chunk
        n_tree_items = read(io, UInt32)
        sfs_n_of_chunks = read(io, UInt32)
        n_tree_chunks = UInt32(ceil((n_tree_items * 0x200) / (chunksize - 0x20)))
        if n_tree_chunks == 1 # file tree do not exceed one chunk in bcf
            seek(io, chunksize * tree_start_chunk + 0x138)
            tree_buffer = IOBuffer(read(io, 0x200 * n_tree_items))
        else
            tree_buffer = IOBuffer()
            tree_pointer = tree_start_chunk
            n_tree_item_chunks = ceil(usable_chunksize / 0x200)
        end

        return version, chunksize, usable_chunksize, tree_start_chunk, n_tree_items, sfs_n_of_chunks, n_tree_chunks, tree_buffer
    end
end

struct SfsItem
    chunk_index_of_pointer_table::UInt32

end
