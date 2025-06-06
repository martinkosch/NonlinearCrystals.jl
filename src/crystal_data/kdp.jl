function create_kdp()
    metadata = Dict(
        :description => "Potassium Dihydrogen Phosphate (KDP)",
        :formula => "KH₂P0₄",
        :point_group => "-42m",
        :lattice_params => (7.448u"Å", 7.448u"Å", 6.977u"Å"),
        :density => 2.3383u"g/cm^3",
        :mohs_hardness => 2.5,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    lambda_min = 0.174u"µm"
    lambda_max = 1.57u"µm"
    temp_ref = 293.15u"K"
    
    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#ghosh1982
    n_o_principal = SellmeierFunction(
        (λ, T) -> begin
            A0 = 1.44896 + 3.185e-5u"K^-1" * T
            B1 = 0.84181 - 1.4114e-4u"K^-1" * T
            C1 = (0.0128 - 2.13e-7u"K^-1" * T) * 1u"µm^2"
            B2 = 0.90793 + 5.75e-7u"K^-1" * T
            C2 = 30.0u"µm^2"

            sqrt(A0 + (B1 * λ^2) / (λ^2 - C1) + (B2 * λ^2) / (λ^2 - C2))
        end,
        (lambda_min, lambda_max)
    )

    n_e_principal = SellmeierFunction(
        (λ, T) -> begin
            A0 = 1.42691 - 1.152e-5u"K^-1" * T
            B1 = 0.72722 - 6.139e-5u"K^-1" * T
            C1 = (0.01213 + 3.104e-7u"K^-1" * T) * 1u"µm^2"
            B2 = 0.22543 - 1.98e-7u"K^-1" * T
            C2 = 30.0u"µm^2"

            sqrt(A0 + (B1 * λ^2) / (λ^2 - C1) + (B2 * λ^2) / (λ^2 - C2))
        end,
        (lambda_min, lambda_max)
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d36=0.39u"pm/V")

    miller_delta = calc_miller_delta(
        d_XYZ_full,
        n_o_principal,
        n_e_principal,
        temp_ref;
        lambda_r1=1064u"nm",
        lambda_r2=1064u"nm",
    )

    KDP = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full;
        miller_delta,
    )
    return KDP
end

export KDP
KDP = create_kdp()