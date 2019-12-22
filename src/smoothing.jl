function buildSavitzkyGolay(coeffs)
        res = Array{Rational{Int}}(undef, 2 * length(coeffs) - 1)
        for i in eachindex(coeffs)
                res[i] = coeffs[i]
                res[end-(i-1)] = coeffs[i]
        end
        return res
end

const SV2 = buildSavitzkyGolay( (17 // 35, 12 // 35, -3 // 35) )
const SV3 = buildSavitzkyGolay( (7 // 21, 6 // 21, 3 // 21, -2 // 21) )
const SV4 = buildSavitzkyGolay( (59 // 231, 54 // 231, 39 // 231, 14 // 231, -21 // 231) )
const SV5 = buildSavitzkyGolay( (89 // 429, 84 // 429, 69 // 429, 44 // 429, 9 // 429, -36 // 429) )
const SV6 = buildSavitzkyGolay( (25 // 143, 24 // 143, 21 // 143, 16 // 143, 9 // 143, 0 // 143, -11 // 143) )
const SV7 = buildSavitzkyGolay( (
        167 // 1105,
        162 // 1105,
        147 // 1105,
        122 // 1105,
        87 // 1105,
        42 // 1105,
        -13 // 1105,
        -78 // 1105,
) )
const SV8 = buildSavitzkyGolay( (
        43 // 323,
        42 // 323,
        39 // 323,
        34 // 323,
        27 // 323,
        18 // 323,
        7 // 323,
        -6 // 323,
        -21 // 323,
) )
const SV9 = buildSavitzkyGolay( (
        269 // 2261,
        264 // 2261,
        249 // 2261,
        224 // 2261,
        189 // 2261,
        144 // 2261,
        89 // 2261,
        24 // 2261,
        -51 // 2261,
        -136 // 2261,
) )
const SV10 = buildSavitzkyGolay( (
        329 // 3059,
        324 // 3059,
        309 // 3059,
        284 // 3059,
        249 // 3059,
        204 // 3059,
        149 // 3059,
        84 // 3059,
        9 // 3059,
        -76 // 3059,
        -171 // 3059,
) )
const SV11 = buildSavitzkyGolay( (
        79 // 805,
        78 // 805,
        75 // 805,
        70 // 805,
        63 // 805,
        54 // 805,
        43 // 805,
        30 // 805,
        15 // 805,
        -2 // 805,
        -21 // 805,
        -42 // 805,
) )
const SV12 = buildSavitzkyGolay( (
        467 // 5175,
        462 // 5175,
        447 // 5175,
        422 // 5175,
        387 // 5175,
        342 // 5175,
        287 // 5175,
        222 // 5175,
        147 // 5175,
        62 // 5175,
        -33 // 5175,
        -138 // 5175,
        -253 // 5175,
) )

function applyfilter(spec::Spectrum, filt, applyLLD = true)
        ty = typeof(promote(spec.counts[1], filt[1])[1])
        sc, fl = counts(spec, ty, applyLLD), convert.(ty, filt)
        res = Array{ty}(undef, length(spec.counts))
        for c = length[filt/2]+1:length(spec)-(length(filt)-length(filt / 2))
                res[c] = sc[c-length[filt/2]:c+length(filt)-length(filt / 2)].filt
        end
        return res
end

filtered(spec::Spectrum, filt, applyLLD = true)::Spectrum =
        Spectrum(spec.energy, applyfilter(spec, filter, applyLLD), copy(spec.properties))
