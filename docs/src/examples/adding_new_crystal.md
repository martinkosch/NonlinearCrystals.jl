# Adding a new nonlinear crystal

To define a new nonlinear crystal in NonlinearCrystals.jl, you should to:

1. Provide metadata such as name, formula, symmetry group, and lattice parameters.
2. Define its principal refractive indices using a `SellmeierFunction` for each crystal axis.
3. Provide the 3×3×3 nonlinear susceptibility tensor `d_XYZ_full` based on individual measured tensor components.
4. Optionally, compute and include a Miller delta correction for frequency scaling.
5. Return the crystal as either a `UnidirectionalCrystal` or `BidirectionalCrystal`.

The definition should be wrapped in a `create_<name>()` function inside a file in `src/crystal_data/`.

## Example: Defining BBO as an unidirectional nonlinear crystal
Below is a complete definition for BBO (Beta-barium borate), a widely used uniaxial nonlinear optical crystal. A similar file can be found in the [crystal_data](https://github.com/martinkosch/NonlinearCrystals.jl/tree/main/src/crystal_data) folder. 

Each crystal has its own file with the wrapping function, e.g.:



```julia
function create_bbo()
    #... All the needed information as shown below should be placed here
    return BBO
end

# Export and instantiate the crystal
export BBO
BBO = create_bbo()
```

Inside the function, follow the steps above.

1. Metadata: Add basic metadata describing material properties as a dict. Especially the point group is needed to calculate the nonlinear tensor matching the crystal symmetry  based on individual components: 
```julia
metadata = Dict(
    :description => "BBO (Beta-barium borate)",
    :formula => "β-BaB₂O₄",
    :point_group => "3m",
    :lattice_params => (12.532u"Å", 12.532u"Å", 12.717u"Å"),
    :density => 3.85u"g/cm^3",
    :mohs_hardness => 4,
)
```

2. Principal refractive indices: Specify the valid spectral range of Sellmeier model and add the refractive data as `SellmeierFunction` with a refractive index function and a reference temperature: 

```julia
    # Specify the valid spectral range and the reference temperature (typically room temperature) of the Sellmeier model
    lambda_min = 0.189u"µm"
    lambda_max = 3.5u"µm"
    temp_ref = 293.15u"K"

    # Define ordinary refractive index as an anonymous function, potentially also incorporating thermo-optic data
    n_o = SellmeierFunction(
        (λ, T) -> sqrt(2.7359 + 0.01878u"µm^2" / (λ^2 - 0.01822u"µm^2") - 0.01354u"µm^-2" * λ^2) - 16.6e-6u"K^-1" * (T - temp_ref),
        (lambda_min, lambda_max);
        temp_ref,
    )

    # Define extraordinary refractive index in a similar way
    n_e = SellmeierFunction(
        (λ, T) -> sqrt(2.3753 + 0.01224u"µm^2" / (λ^2 - 0.01667u"µm^2") - 0.01516u"µm^-2" * λ^2) - 9.3e-6u"K^-1" * (T - temp_ref),
        (lambda_min, lambda_max);
        temp_ref,
    )
```

3. The nonlinear susceptibility tensor: 
The function `d_XYZ_full` lets you compute the full nonlinear tensor in the dielectric frame. Rotations into this frame can be specified as a rotation matrix if needed. See the [file for BIBO](https://github.com/martinkosch/NonlinearCrystals.jl/blob/main/src/crystal_data/bibo.jl) as an example.
In the present case, `d15` and `d22` are the two non-zero coefficients for BBO's `3m` point group; all other components are set automatically based on the crystal's point group: 
```julia
    d_XYZ_full = calc_d_XYZ_full(
        metadata[:point_group]; 
        d22 = -2.2u"pm/V", 
        d15 = 0.08u"pm/V"
    )
```
The full or contracted tensors can in the end be shown using `BBO.d_XYZ_ref_full` or `BBO.d_XYZ_ref`, respectively. 

4. Miller delta: Compute Miller scaling data for corrections of `d_eff` if the wavelengths during the measurements of the specified tensor components are known. For example, if the specified tensor components `d15` and `d22` are known to be measured using the second-harmonic generation of 1064 nm light, calculate the miller delta tensor using: 
```julia
    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_o, 
        n_e, 
        temp_ref; 
        lambda_r1 = 1064u"nm", 
        lambda_r2 = 1064u"nm"
    )
```

5. Returning the crystal object: Finally, construct and return the actual crystal object:
```julia
    BBO = UnidirectionalCrystal(
        metadata,
        n_o,
        n_e,
        d_XYZ_full;
        miller_delta,
    )
    return BBO
```

It is helpful to add the references to the publications of the used data as comments. 
Make also sure that the specified data is working as expected by adding a corresponding test file for each new nonlinear crystal. They should be placed in the [crystal data test folder](https://github.com/martinkosch/NonlinearCrystals.jl/tree/main/test/crystal_data_tests). 