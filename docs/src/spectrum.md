# ![](NeXL_sm.png)Spectrum
## Working with Spectrum objects

```@meta
CurrentModule = NeXLSpectrum
```

The `Spectrum` type represents a single spectrum with associated properties.
The channel data is indexed like a `Vector`.  The property data is indexed
using `Symbol` objects.  A set of spectrum properties are defined in the
library and the user can create additional ones.

In addition to the `Vector`-like properties, the `Spectrum` associates
energy bins with each channel.  The energy bins are continuous, non-overlapping
and monotonic.  Functions like `energy(...)` maps channel index to energies
and `channel(...)` to make energies to channel index.  

```@docs
Spectrum
```

# Acccessing Properties
Most properties are accessed using `spec[:Symbol]` notation.  Some combined
or special properties have special methods.
```@docs
dose
elms(::Spectrum)
```

# Displaying Spectrum Data
```@docs
asa
Base.keys
Base.haskey
Gadfly.plot
```

# Energy Scale Functions
These functions handle mapping channel index to energies and vice versa.

```@docs
NeXLCore.energy
channel
rangeofenergies
channelwidth
energyscale
```

# Extracting the Counts Data
The channel data in a `Spectrum` object `spec` can be get/set using standard
`Array` indexing techniques.

    > spec[123]
    > spec[123:345]
    > spec[1234.0] # Channel containing energy 1234.0 eV
    > spec[123] = 99

Alternatively, the `counts(...)` method can be used.

```@docs
counts
```

`counts(...)` can apply the :LLD property to zero counts below a specified channel.
```@docs
lld
```

# Defining Detectors
Build a detector to match the data in a `Spectrum`.
```@docs
simpleEDS
matching
```

# Statistical Sub-Sampling
```@docs
NeXLSpectrum.subdivide
NeXLSpectrum.subsample

```

# Processing Channel Data
```@docs
kratio
integrate
estimatebackground
modelBackground
extractcharacteristic
peak
background
peaktobackground
estkratio
```

# Miscellaneous Functions
```@docs
details
commonproperties
maxspectrum
Base.findmax
```
