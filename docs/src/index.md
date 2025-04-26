```@meta
CurrentModule = NonlinearCrystals
```

# NonlinearCrystals.jl

NonlinearCrystals.jl is a Julia package for analyzing and simulating nonlinear optical processes in birefringent crystals. It provides tools for evaluating refractive data, finding phase-matches, and generating interactive plots of critical and noncritical phase-matches.

## Installation

This package is registered in the general Julia registry. To install the package, use Julia's package manager:

```julia
pkg> add NonlinearCrystals
```

## Current Features

- **Refraction modeling** for uniaxial and biaxial nonlinear crystals  
- Computation of **group velocities**, **walkoff angles**, and **dispersion terms** 
- Calculation of **effective nonlinearity** $d_\text{eff}$ with or without **Miller scaling**
- Determination of critical and noncritical **collinear phase-matching** conditions
- Strict use of **units** based on [Unitful.jl](https://github.com/PainterQubits/Unitful.jl/tree/master)
- **Visualization tools** for:
  - Wavelength and temperature-dependent refractive indices
  - Miller scaling of nonlinear coefficients
  - Various Δk plots for critical phase-matches 
  - Noncritical phase-matching lines and SHG points