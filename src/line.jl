# Line drawing algorithm

function drawline(func::Function, ci1::CartesianIndex{2}, ci2::CartesianIndex{2}, eachstep::Bool=false)
  r, c = ci1.I
  dr, dc = -abs(ci2.I[1] - ci1.I[1]), abs(ci2.I[2] - ci1.I[2])
  sr, sc = ci1.I[1] < ci2.I[1] ? 1 : -1, ci1.I[2] < ci2.I[2] ? 1 : -1
  err = dc + dr
  func( ci1.I )
  while !((c == ci2.I[2]) && (r == ci2.I[1]))
    e2 = 2*err
    if e2 >= dr
      err += dr
      c += sc
      if eachstep
        func( ( r, c ) )
      end
    end
    if e2 <= dc
      err += dc
      r += sr
      if eachstep
        func( ( r, c ) )
      end
    end
    if !eachstep
        func( ( r, c ) )
    end
  end
end
