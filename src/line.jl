# Line drawing algorithm

"""
    drawline(func::Function, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}, eachstep::Bool=false)

At each step along the line from `ci1` to `ci2` call `func` with a single argument, the row, column coordinates
of the pixel as a `Tuple{Int, Int}`.

If `eachstep=true` then `func` is called once when the row changes and once when the column changes.  When
`eachstep=false` then `func` is called less frequently and both row and col may change between calls.
"""
function drawline(
    func::Function,
    ci1::CartesianIndex{2},
    ci2::CartesianIndex{2},
    eachstep::Bool = false,
)
    r, c = ci1.I
    dr, dc = -abs(ci2.I[1] - ci1.I[1]), abs(ci2.I[2] - ci1.I[2])
    sr, sc = ci1.I[1] < ci2.I[1] ? 1 : -1, ci1.I[2] < ci2.I[2] ? 1 : -1
    err = dc + dr
    func(ci1.I)
    while !((c == ci2.I[2]) && (r == ci2.I[1]))
        e2 = 2 * err
        if e2 >= dr
            err += dr
            c += sc
            if eachstep
                func((r, c))
            end
        end
        if e2 <= dc
            err += dc
            r += sr
            if eachstep
                func((r, c))
            end
        end
        if !eachstep
            func((r, c))
        end
    end
end

"""
    drawray(
        func::Function,  # returns Bool 
        ci1::CartesianIndex{2}, # starting point
        ci2::CartesianIndex{2}, # towards this point
        bounds::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}},  # stay within
        eachstep::Bool=false
    )
    function drawray(
        func::Function,
        ci1::CartesianIndex{2},
        ci2::CartesianIndex{2},
        eachstep::Bool = false,
    )

    Draw a ray starting at `ci1` towards `ci2` and call `func` with a single argument, the row, column coordinates
of the pixel as a `Tuple{Int, Int}`.  `func` must return a `Bool` - true to continue drawing the ray and false to
stop drawing the ray.  The first version checks the bounds each step and the second doesn't.

If `eachstep=true` then `func` is called once when the row changes and once when the column changes.  When
`eachstep=false` then `func` is called less frequently and both row and col may change between calls.
"""

function drawray(
    func::Function,
    ci1::CartesianIndex{2},
    ci2::CartesianIndex{2},
    bounds::CartesianIndices{2, Tuple{Vararg{OrdinalRange{Int64, Int64}, 2}}},
    eachstep::Bool = false,
)
    r, c = ci1.I
    dr, dc = -abs(ci2.I[1] - ci1.I[1]), abs(ci2.I[2] - ci1.I[2])
    sr, sc = ci1.I[1] < ci2.I[1] ? 1 : -1, ci1.I[2] < ci2.I[2] ? 1 : -1
    err = dc + dr
    inb(r,c) = ((r in bounds.indices[1]) && (c in bounds.indices[2]))
    if !(inb(r,c) && func((r,c)))
        return nothing
    end
    while true  
        e2 = 2 * err
        if e2 >= dr
            err += dr
            c += sc
            if eachstep && !(inb(r,c) && func((r, c)))
                return nothing
            end
        end
        if e2 <= dc
            err += dc
            r += sr
            if eachstep && !(inb(r,c) && func((r, c)))
                return nothing
            end
        end
        if (!eachstep) && !(inb(r,c) && func((r, c)))
            return nothing
        end

    end
end

function drawray(
    func::Function,
    ci1::CartesianIndex{2},
    ci2::CartesianIndex{2},
    eachstep::Bool = false,
)
    r, c = ci1.I
    dr, dc = -abs(ci2.I[1] - ci1.I[1]), abs(ci2.I[2] - ci1.I[2])
    sr, sc = ci1.I[1] < ci2.I[1] ? 1 : -1, ci1.I[2] < ci2.I[2] ? 1 : -1
    err = dc + dr
    if !func((r,c))
        return nothing
    end
    while true  
        e2 = 2 * err
        if e2 >= dr
            err += dr
            c += sc
            if eachstep && !func((r, c))
                return nothing
            end
        end
        if e2 <= dc
            err += dc
            r += sr
            if eachstep && !func((r, c))
                return nothing
            end
        end
        if (!eachstep) && !func((r, c))
            return nothing
        end

    end
end