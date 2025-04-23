# NonlinearCrystals.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinkosch.github.io/NonlinearCrystals.jl/dev/)
[![Build Status](https://github.com/martinkosch/NonlinearCrystals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/martinkosch/NonlinearCrystals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/martinkosch/NonlinearCrystals.jl/graph/badge.svg?token=AVIIA69G5F)](https://codecov.io/gh/martinkosch/NonlinearCrystals.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

NonlinearCrystals.jl is a Julia package for analyzing and simulating nonlinear optical processes in birefringent crystals. It provides tools for evaluating refractive data, finding phasematches, and generating plots of critical and noncritical phase matches.

## Installation

This package is not yet registered in the general Julia registry. To install the package, clone it using Julia's package manager:

```julia
]dev https://github.com/martinkosch/NonlinearCrystals.jl
```

This will clone the package locally with its latest development version.

## Current features

- **Refraction modeling** for uniaxial and biaxial nonlinear crystals  
- Computation of **group velocities**, **walkoff angles**, and **dispersion terms** 
- Calculation of **effective nonlinearity** $d_\text{eff}$ with or without **Miller scaling**
- Determination of critical and noncritical **collinear phase-matching** conditions
- Strict use of **units** based on [Unitful.jl](https://github.com/PainterQubits/Unitful.jl/tree/master)
- **Interactive visualization** for:
  - Wavelength and temperature-dependent refractive indices
  - Miller scaling of nonlinear coefficients
  - Various Δk plots for critical phasematches 
  - Noncritical phasematching lines and SHG points

## Example usage

You can use the package to inspect the nonlinear and refractive properties of various crystals, find phase-matching directions, and generate interactive plots.

For example, to find and analyze a critical phasematch for third-harmonic generation of 1064 nm light in LBO, one can first plot all possible critical phasematches:

```julia
using NonlinearCrystals
using Unitful

# Define phasematch parameters
lambda_r1=1064u"nm"
lambda_r2=532u"nm"
temp = 293.15u"K"
hi_or_lo_rrb=[(:hi, :hi, :lo), (:lo, :lo, :hi)]

# Display and plot all possible phasematches with the specified combination of polarizations
plot_critical_pms(LBO; hi_or_lo_rrb, lambda_r1, lambda_r2, temp)
```

![Screenshot of all critical phasematches returned by the plot_critical_pms function.](https://github.com/martinkosch/NonlinearCrystals.jl/blob/main/docs/src/lbo_all_pms.png)

It is also possible to visualize the phasematches for a given combination of polarizations on a unit sphere representing possible incidence directions. 
This plot shows details on the selected critical phasematch as a label on mouseover. This functionality might take a few seconds for precompilation. Clicking on a point on the phasematch contour prints the phasematch to the terminal and appends it to the exported global variable `selected_pms` for later use.  
```julia
# Visualize critical phase matches for (:hi, :hi, :lo) on a unit sphere
plot_delta_k_map(hi_or_lo_rrb[1], LBO; lambda_r1, lambda_r2, temp, plot_type=:sphere)
```
![Screenshot of the phase mismatch Δk over all incidence directions returned by the plot_delta_k_map function.](https://github.com/martinkosch/NonlinearCrystals.jl/blob/main/docs/src/lbo_example.png)

To select a certain phasematch, various search methods are available. For example, to search for a critical phasematch close to a certain incidence direction, one can use the following function: 
```julia
# Find a specific phasematch
theta = 90u"°"
phi = 30u"°"
julia> pm = find_nearest_pm_along_theta_phi(theta, phi, hi_or_lo_rrb[1], LBO; lambda_r1, lambda_r2, temp)

Crystal:                      LBO (Lithium Triborate)
k angles:                     θ: 90.00°, ϕ: 37.22°
k direction:                  [0.796, 0.605, 0.0]      
Temperature:                  293.15 K (20.00 °C)
────────────────────────────────────────────────────────────────────────────────────────────────────────
Wavelength (nm):              1064.0                    532.0                     354.667                  
Refractive index type:        hi                        hi                        lo                       
Type I PM in XY plane:        o                         o                         e                        
Phase velocity / c₀:          1.605                     1.622                     1.616                    
Group velocity / c₀:          1.626                     1.656                     1.690                    
Walkoff angle (mrad):         0.000                     0.000                     18.070                   
S direction:                  [0.796, 0.605, 0.0]       [0.796, 0.605, 0.0]       [0.807, 0.59, 0.0]       
E direction:                  [-0.0, -0.0, 1.0]         [-0.0, -0.0, 1.0]         [0.59, -0.807, 0.0]      
D direction:                  [-0.0, -0.0, 1.0]         [-0.0, -0.0, 1.0]         [0.605, -0.796, 0.0]     
GDD (fs²/mm):                 16.592                    89.684                    153.153                  
TOD (fs³/mm):                 70.119                    35.846                    46.928                   
ω BW × L (GHz·cm):            872.759                   464.602                   993.454                  
T BW × L (K·cm):              20.621                   
θ BW × L (mrad·cm):           Inf                      
ϕ BW × L (mrad·cm):           1.202                    
d_eff (pm/V):                 -0.704 (w/o Miller scaling: -0.686)
S₀ × L² (W):                  1.620e+08 (w/o Miller scaling: 1.700e+08)
────────────────────────────────────────────────────────────────────────────────────────────────────────
```
More examples are available in the [examples](https://martinkosch.github.io/NonlinearCrystals.jl/dev/) section.

## Roadmap

Development is actively ongoing. Planned features include:

- Improved **interactive visualizations**
- **Quasi-phasematching (QPM)** support  
- 1D and 2D static and transient **simulations** of nonlinear three-wave processes
- **Non-collinear** phasematches
- More **crystal data**
