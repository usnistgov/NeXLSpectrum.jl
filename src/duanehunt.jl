# Duane-Hunt related functions

_duane_hunt_func(es, p) = map(ee-> p[1]*bremsstrahlung(Small1987, ee, p[2], n"H"), es)

function _duane_hunt_impl(spec::Spectrum)
    if channel(spec[:BeamEnergy], spec) < length(spec)
        cdata = counts(spec)
        este0, estCh = let # Sometimes the D-H can be much below `spec[:BeamEnergy]`
            # this handles pulse pile-up above the D-H
            lim = min(10.0, max(1.0, 1.0e-5 * maximum(cdata)))
            len = min(10, max(5, length(spec) ÷ 400))
            chlast = findlast(i -> mean(@view cdata[i-len:i]) > lim, len+1:length(cdata))
            energy(chlast, spec), chlast
        end
        xdata, ydata = let
            # Range of channels above/below `este0`
            chs = max(1,(9*estCh)÷10):min((21*estCh)÷20,length(spec))
            energyscale(spec, chs), cdata[chs]
        end
        esti0 = sum(ydata) / sum(_duane_hunt_func(xdata, (1.0, este0)))
        e0 = spec[:BeamEnergy]
        # Apply constraint functions to keep DH between 0.1 E0 and 1.1 E0.
        transform(v) = asin(2*(v-0.1*e0)/(1.1*e0-0.1*e0)-1)
        inv_transform(tv) = 0.1*e0+(sin(tv)+1)*(1.1*e0-0.1*e0)/2
        # Fit Small's continuum model to the data (ignore detector efficiency)
        model(xs, p) = _duane_hunt_func(xs, (p[1], inv_transform(p[2])))
        res = curve_fit(model, xdata, ydata, [esti0, transform(este0)])
        return (res.param[1], inv_transform(res.param[2]))
    else
        error("Unable to estimate the Duane-Hunt on $(spec[:Name]) because the counts data does not extend to $(spec[:BeamEnergy]) keV.")
    end
end

"""
    duane_hunt(spec::Spectrum)

Estimates the Duane-Hunt limit (the energy at which the continuum goes to effectively zero.)
"""
duane_hunt(spec::Spectrum)  = _duane_hunt_impl(spec)[2]