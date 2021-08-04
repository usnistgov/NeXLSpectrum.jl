## Quantifying K412 using NeXLSpectrum and NeXLMatrixCorrection

This document demonstrates the low-level API for filter fitting spectra.  It is more flexible than the higher-
level API but much more complex.  In most situations, the high level API discussed [here](K412refs.md) is more appropriate.

Use the NeXLSpectrum to load, plot, fit and report the quantification of a set of K412 spectra.

Loading `NeXLSpectrum` also automatically makes `NeXLCore` and `NeXLUncertainties` available.

Loading the `Gadfly` library adds plotting support to `NeXLSpectrum`.

```julia
using NeXLSpectrum              # Provides spectrum reading and fitting tools
using NeXLMatrixCorrection      # Provides `quant` to convert k-ratios to mass fraction.
using Gadfly                    # Plotting
using DataFrames, Latexify      # Tables
```



#### Read in the Spectra
```julia
path = joinpath(@__DIR__,"K412 spectra")
# Load a single spectrum
fe = loadspectrum(joinpath(path, "Fe std.msa"))
# Create a detector model to match it
det = matching(fe, 132.0, 10)
# Now load all the spectra using this detector
unks = (i->loadspectrum(joinpath(path, "III-E K412[$i][4].msa"),det)).(0:4)
al2o3 = loadspectrum(joinpath(path, "Al2O3 std.msa"),det)
caf2 = loadspectrum(joinpath(path, "CaF2 std.msa"),det)
fe = loadspectrum(joinpath(path, "Fe std.msa"),det)
mgo = loadspectrum(joinpath(path, "MgO std.msa"),det)
sio2 = loadspectrum(joinpath(path, "SiO2 std.msa"),det)
# Add carbon coating
foreach(s->s[:Coating]=Film(pure(n"C"), 30.0e-7), unks)
foreach(s->s[:Coating]=Film(pure(n"C"), 10.0e-7), (al2o3, caf2, mgo, sio2));
```




Table: The spectra

|               Name | BeamEnergy | ProbeCurrent | LiveTime | RealTime |           Coating |  Integral | Material |
| ------------------:| ----------:| ------------:| --------:| --------:| -----------------:| ---------:| --------:|
| III-E K412[0][all] |      2e+04 |        1.114 |    235.5 |    286.3 | 30.0 nm of Pure C | 8.079e+06 |     K412 |
| III-E K412[1][all] |      2e+04 |        1.114 |    235.4 |    286.2 | 30.0 nm of Pure C | 8.077e+06 |     K412 |
| III-E K412[2][all] |      2e+04 |        1.112 |    235.5 |    286.3 | 30.0 nm of Pure C | 8.084e+06 |     K412 |
| III-E K412[3][all] |      2e+04 |         1.11 |    235.4 |    286.3 | 30.0 nm of Pure C | 8.087e+06 |     K412 |
| III-E K412[4][all] |      2e+04 |         1.11 |    235.4 |    286.2 | 30.0 nm of Pure C | 8.081e+06 |     K412 |
|          Al2O3 std |      2e+04 |         1.11 |     1172 |     1491 | 10.0 nm of Pure C | 4.974e+07 |    Al2O3 |
|           CaF2 std |      2e+04 |         1.11 |     1176 |     1456 | 10.0 nm of Pure C | 4.406e+07 |     CaF2 |
|             Fe std |      2e+04 |         1.11 |     1171 |     1529 |           nothing | 5.445e+07 |       Fe |
|            MgO std |      2e+04 |        1.106 |     1176 |     1496 | 10.0 nm of Pure C | 4.985e+07 |      MgO |
|           SiO2 std |      2e+04 |         1.11 |     1173 |     1470 | 10.0 nm of Pure C | 4.665e+07 |     SiO2 |



Notice that the spectra all have 1) live-time (`:LiveTime`); 2) probe-current (`:ProbeCurrent`); 3) take-off angle
(`:TakeOffAngle`); 4) beam energy (`:BeamEnergy`); and detector (`:Detector`) properties defined.  These properties
are necessary for extracting the k-ratios and estimating the composition.
```julia
sio2[:LiveTime], sio2[:ProbeCurrent], sio2[:TakeOffAngle], sio2[:BeamEnergy], sio2[:Detector]
```

```
(1173.1648, 1.10989, 0.6108652381980153, 20000.0, BasicEDS[4096 chs, E[ch] 
= 1.63032 + 9.99856⋅ch, 132.0 eV @ Mn K-L3, 10 ch LLD, [Be,Sc,Ba,Pu]])
```




#### The Unknowns
```julia
display(plot(unks..., klms=[n"O",n"Mg",n"Al",n"Si",n"Ca",n"Fe"], xmax=8.0e3))
```

![](figures/K412fit_5_1.svg)


#### The Reference Spectra
Build a convenient structure so it is easy to appreciate the necessary information and to splat it into
`filteredReference`.
```julia
refs = (
  # spectrum, element, composition
  ( al2o3, n"Al", mat"Al2O3" ), #
  ( mgo,   n"Mg", mat"MgO" ),   #
  ( fe,    n"Fe", mat"Fe" ),    #
  ( sio2,  n"Si", mat"SiO2" ),  #
  ( sio2,  n"O",  mat"SiO2" ),  #
  ( caf2,  n"Ca", mat"CaF2" ), )
display(plot(al2o3, caf2, fe, mgo, sio2, klms=collect( ref[2] for ref in refs), xmax=8.0e3))
```

![](figures/K412fit_6_1.svg)



#### Pre-filter the Reference Spectra
```julia
# Build a top-hat filter
filt = buildfilter(NeXLSpectrum.GaussianFilter,det)
# Filter all the reference spectra
frs = mapreduce(ref->filterreference(filt, ref..., withEsc=true), append!, refs)
# frs is now a FilteredReference[] used to fit the unknowns.
```

```
12-element Vector{FilteredReference}:
 Reference[k[Al K-L3 + 1 other, Al2O3]]
 Reference[k[Mg K-L3 + 1 other, MgO]]
 Reference[k[Fe L3-M5 + 11 others, Fe]]
 Reference[k[Fe K-L3 + 1 other, Fe]]
 Reference[k[Fe K-M3 + 3 others, Fe]]
 Reference[Ecs[Fe K-L3 + 1 other]]
 Reference[Ecs[Fe K-M3 + 3 others]]
 Reference[k[Si K-L3 + 2 others, SiO2]]
 Reference[k[O K-L3 + 1 other, SiO2]]
 Reference[k[Ca K-L3 + 3 others, CaF2]]
 Reference[Ecs[Ca K-L3 + 1 other]]
 Reference[Ecs[Ca K-M3 + 1 other]]
```





#### Fit the Pre-Filtered References to the Unknowns
```julia
res= [ fit_spectrum(unk,filt,frs,false) for unk in unks ]
```

```
5-element Vector{FilterFitResult}:
 III-E K412[0][all]
 III-E K412[1][all]
 III-E K412[2][all]
 III-E K412[3][all]
 III-E K412[4][all]
```




|            Spectra | k[O K-L3 + 1 other, SiO2] | k[Fe L3-M5 + 11 others, Fe] | k[Mg K-L3 + 1 other, MgO] | k[Al K-L3 + 1 other, Al2O3] | k[Si K-L3 + 2 others, SiO2] | k[Ca K-L3 + 3 others, CaF2] | k[Fe K-L3 + 1 other, Fe] | k[Fe K-M3 + 3 others, Fe] |
| ------------------:| -------------------------:| ---------------------------:| -------------------------:| ---------------------------:| ---------------------------:| ---------------------------:| ------------------------:| -------------------------:|
| III-E K412[0][all] |                     0.652 |                     0.04245 |                    0.1475 |                     0.06704 |                       0.351 |                      0.1922 |                  0.06684 |                   0.06675 |
| III-E K412[1][all] |                     0.654 |                     0.04208 |                    0.1474 |                     0.06677 |                      0.3501 |                      0.1916 |                  0.06706 |                   0.06733 |
| III-E K412[2][all] |                    0.6545 |                     0.04245 |                    0.1478 |                     0.06713 |                      0.3513 |                      0.1922 |                  0.06686 |                   0.06691 |
| III-E K412[3][all] |                    0.6589 |                     0.04199 |                     0.148 |                      0.0672 |                      0.3521 |                      0.1926 |                  0.06682 |                   0.06768 |
| III-E K412[4][all] |                    0.6573 |                     0.04136 |                    0.1481 |                     0.06731 |                       0.352 |                      0.1922 |                  0.06692 |                   0.06648 |




Let's take a look at a residual spectrum by plotting one of the `FilterFitResult` objects.
```julia
plot(res[1])
```

![](figures/K412fit_10_1.svg)



#### Quantify the k-ratios by Matrix Correction
```julia
quant = quantify.(res)
```

```
5-element Vector{IterationResult}:
 Converged to III-E K412[0][all][Al=0.0487,Ca=0.1089,Fe=0.0807,Mg=0.1175,Si
=0.2084,O=0.4795] in 9 steps.
 Converged to III-E K412[1][all][Al=0.0485,Ca=0.1085,Fe=0.0809,Mg=0.1175,Si
=0.2079,O=0.4800] in 10 steps.
 Converged to III-E K412[2][all][Al=0.0488,Ca=0.1089,Fe=0.0807,Mg=0.1178,Si
=0.2086,O=0.4809] in 9 steps.
 Converged to III-E K412[3][all][Al=0.0488,Ca=0.1091,Fe=0.0807,Mg=0.1179,Si
=0.2090,O=0.4835] in 10 steps.
 Converged to III-E K412[4][all][Al=0.0489,Ca=0.1089,Fe=0.0808,Mg=0.1180,Si
=0.2090,O=0.4824] in 10 steps.
```




|           Material |      O |     Mg |      Al |     Si |     Ca |      Fe | Total |
| ------------------:| ------:| ------:| -------:| ------:| ------:| -------:| -----:|
| III-E K412[0][all] | 0.4795 | 0.1175 |  0.0487 | 0.2084 | 0.1089 | 0.08067 | 1.044 |
| III-E K412[1][all] |   0.48 | 0.1175 | 0.04852 | 0.2079 | 0.1085 | 0.08093 | 1.043 |
| III-E K412[2][all] | 0.4809 | 0.1178 | 0.04877 | 0.2086 | 0.1089 | 0.08069 | 1.046 |
| III-E K412[3][all] | 0.4835 | 0.1179 | 0.04882 |  0.209 | 0.1091 | 0.08065 | 1.049 |
| III-E K412[4][all] | 0.4824 |  0.118 |  0.0489 |  0.209 | 0.1089 | 0.08077 | 1.048 |




Finally plot the results as mass fractions.
```julia
plot(quant, known=unks[1][:Composition])
```

![](figures/K412fit_13_1.svg)



Plot the difference from the SRM value.
```julia
plot(quant, known=unks[1][:Composition], delta=true)
```

![](figures/K412fit_14_1.svg)



Plot the difference from the mean value for each element.
```julia
plot(quant, delta=true)
```

![](figures/K412fit_15_1.svg)
