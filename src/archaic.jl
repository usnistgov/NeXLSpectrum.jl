"""
    function labeledextents(
        elm::Element,  # All CharXRay for this element
        det::Detector,
        ampl::Float64,
        maxE::Float64=energy(det.channelcount+1, det) # full detector range
    )::Vector{Tuple{Vector{CharXRay},UnitRange{Int}}}

Creates a vector containing pairs containing a vector of CharXRay and an interval. The interval represents a
contiguous interval over which all the X-rays in the interval are sufficiently close in energy that they will
interfere with each other on the specified detector.
"""
labeledextents(elm::Element, det::Detector, ampl::Float64, maxE::Float64=1.0e6) =
    labeledextents(visible(characteristic(elm, alltransitions, ampl, maxE), det), det, ampl)

"""
    function labeledextents(
        elm::Vector{Element},  # All CharXRay for these elements
        det::Detector,
        ampl::Float64,
        maxE::Float64=energy(det.channelcount+1, det)
    )::Vector{Tuple{Vector{CharXRay},UnitRange{Int}}}

Creates a vector containing pairs containing a vector of CharXRay and an interval. The interval represents a
contiguous interval over which all the X-rays in the interval are sufficiently close in energy that they will
interfere with each other on the specified detector.
"""
labeledextents(elms::Vector{Element}, det::Detector, ampl::Float64, maxE::Float64=energy(det.channelcount+1, det)) =
    labeledextents(mapreduce(elm->characteristic(elm,alltransitions),append!,elms),det,ampl)

extents(elm::Element, det::Detector, ampl::Float64)::Vector{UnitRange{Int}} =
    extents(visible(characteristic(elm,alltransitions),det),det,ampl)

extents(elms::Vector{Element}, det::Detector, ampl::Float64)::Vector{UnitRange{Int}} =
    extents(mapreduce(elm->characteristic(elm,alltransitions),append!,elms),det,ampl)
