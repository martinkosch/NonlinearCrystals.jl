function create_gs()
    metadata = Dict(
        :description => "GaSe (Gallium(II) selenide)",
        :formula => "GaSe",
        :point_group => "-6m2",
        :lattice_params => (3.755u"Å", 3.755u"Å", 15.94u"Å"),
        :density => 5.03u"g/cm^3",
        :mohs_hardness => 3,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#vodopyanov1995
    lambda_min = 0.65u"µm" 
    lambda_max = 18u"µm" 
    temp_ref = 293u"K"
    
    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#vodopyanov1995
    n_o_principal = SellmeierFunction(
        (λ, T) -> begin
            n_sq = 7.443 +
                0.4050u"µm^2" / λ^2 +
                0.0186u"µm^4" / λ^4 +
                0.0061u"µm^6" / λ^6 +
                3.1485 * λ^2 / (λ^2 - 2194u"µm^2")
                return sqrt(n_sq)  # Add temperature dependence here if known
        end,
        (lambda_min, lambda_max);
        temp_ref,
        )
        
    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#vodopyanov1995
    n_e_principal = SellmeierFunction(
        (λ, T) -> begin
            n_sq = 5.76 +
                0.3879u"µm^2" / λ^2 -
                0.2288u"µm^4" / λ^4 +
                0.1223u"µm^6" / λ^6 +
                1.8550 * λ^2 / (λ^2 - 1780u"µm^2")
            return sqrt(n_sq)  # Add temperature dependence here if known
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#abdullaev1972
    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d22=54u"pm/V")

    miller_delta = calc_miller_delta(
        d_XYZ_full,
        n_o_principal,
        n_e_principal,
        temp_ref;
        lambda_r1=10.6u"µm",
        lambda_r2=10.6u"µm",
    )

    GS = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full;
        miller_delta,
    )
    return GS
end

export GS
const GS = create_gs()