# Inspecting and plotting the refractive index

NonlinearCrystals.jl provides access to refractive index data for a wide range of birefringent nonlinear crystals. Built-in crystals like `LBO`, `BBO`, or `KTP_F` are readily available as objects. Each crystal is defined in its own file under the [crystal_data](https://github.com/martinkosch/NonlinearCrystals.jl/tree/main/src/crystal_data) directory.

Depending on the nature of the crystal—**unidirectional** or **bidirectional**—its birefringence is described via two or three principal refractive indices that depend on both wavelength and temperature. These are provided as `SellmeierFunction` objects and can be accessed as fields of the crystal:

- For **unidirectional** crystals (e.g., `BBO`):
  ```julia
  BBO.n_o_principal   # ordinary axis
  BBO.n_e_principal   # extraordinary axis
  ```

- For **bidirectional** crystals (e.g., `KTP_F`):
  ```julia
  KTP_F.n_X_principal
  KTP_F.n_Y_principal
  KTP_F.n_Z_principal
  ```

Here's how to evaluate basic refractive and dispersion properties for a principal axis:

```julia
julia> using NonlinearCrystals, Unitful

julia> n = refractive_index(LBO.n_X_principal, 1064u"nm", 300u"K")   # Refractive index for polarization along the principal X axis at 1064 nm and 300 K
1.564762196788999

julia> ng = group_index(LBO.n_X_principal, 1064u"nm", 300u"K")   # Group index
1.58131698872782

julia> β₂ = β2(LBO.n_X_principal, 1064u"nm", 300u"K")   # Group delay dispersion (GDD)
1.7872825157935612e-26 s^2 m^-1
```

All inputs and outputs of NonlinearCrystals.jl have physical units based on Unitful.jl. It is easy to convert those if needed:

```julia
julia> β₂ |> u"fs^2/mm"
17.872825157935612 fs^2 mm^-1
``` 

To analyze thermal effects, you can also compute the thermo-optic coefficient:
```julia
julia> dn_dT = dn_dtemp(LBO.n_X_principal, 1064u"nm", 300u"K")  # ∂n/∂T at 1064 nm and 300 K
-2.3793331118400023e-6 K^-1
```

## Querying refractive indices

These principal indices apply when the electric field is aligned with the corresponding dielectric axis. For arbitrary propagation directions or polarizations, you can compute direction-dependent refractive indices, specified for polar angles θ and ϕ:
```julia
julia> RefractionDataHiLo(45u"°", 10u"°", LNB_C, 532u"nm"; temp= 349u"K")

Crystal:                  LNB (Lithium Niobate, congruent melt with mole ratio Li/Nb 0.946)
Wavelength (nm):          532.0
k angles:                 θ: 45.00°, ϕ: 10.00°
k direction:              [0.696, 0.123, 0.707]    
Temperature:              349.00 K (75.85 °C)
───────────────────────────────────────────────────────────────────────────
Refractive index type:    hi (o)                    lo (e)                   
Phase velocity / c₀:      2.324                     2.279                    
Group velocity / c₀:      2.586                     2.523                    
Walkoff angle (mrad):     0.000                     37.905                   
S direction:              [0.696, 0.123, 0.707]     [0.669, 0.118, 0.733]    
E direction:             ±[-0.174, 0.985, -0.0]    ±[-0.722, -0.127, 0.68]  
D direction:             ±[-0.174, 0.985, -0.0]    ±[-0.696, -0.123, 0.707] 
GDD (fs²/mm):             896.823                   825.030                  
TOD (fs³/mm):             544.009                   489.841
```

## Plotting refractive index vs. wavelength
You can visualize the dispersion behavior using:
```julia
plot_refractiveindex(LBO) 
```

To see the temperature dependence:
```julia
plot_refractiveindex(LBO; temp=[293.15u"K", 413.15u"K", 533.15u"K"])
```
This will generate interactive plots of refractive index as a function of wavelength at different temperatures.