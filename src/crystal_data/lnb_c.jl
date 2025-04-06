function create_lnb_c()
    metadata = Dict(
        :description => "LNB (Lithium Niobate, congruent melt with mole ratio Li/Nb 0.946)",
        :formula => "LiNbO₃",
        :point_group => "3m",
        :lattice_params => (5.148u"Å", 5.148u"Å", 3.863u"Å"),
        :density => 4.628u"g/cm^3",
        :mohs_hardness => 5,
    )

    # Handbook of Nonlinear Crystals
    lambda_min = 0.4u"µm"
    lambda_max = 5.5u"µm"
    temp_ref = 293.15u"K"

    n_o_principal = SellmeierFunction(
        (λ, T) -> sqrt(
            4.9048 +
            2.1429e-8u"K^-2" * (T^2 - 88506.25u"K^2") +
            (0.11775u"µm^2" + 2.2314e-8u"µm^2 * K^-2" * (T^2 - 88506.25u"K^2")) /
            (λ^2 - (0.21802u"µm" - 2.9671e-8u"µm * K^-2" * (T^2 - 88506.25u"K^2"))^2) -
            0.027153u"µm^-2" * λ^2
        ),
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_e_principal = SellmeierFunction(
        (λ, T) -> sqrt(
            4.5820 +
            2.2971e-7u"K^-2" * (T^2 - 88506.25u"K^2") +
            (0.09921u"µm^2" + 5.2716e-8u"µm^2 * K^-2" * (T^2 - 88506.25u"K^2")) /
            (λ^2 - (0.21090u"µm" - 4.9143e-8u"µm * K^-2" * (T^2 - 88506.25u"K^2"))^2) -
            0.021940u"µm^-2" * λ^2
        ),
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d22=2.10u"pm/V", d31=-4.35u"pm/V", d33=-27.2u"pm/V") # Measured at 1.06 µm

    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_o_principal,
        n_e_principal,
        temp_ref;  # TODO: This is a test/guess!
        lambda_r1=800u"nm", # TODO: This is a test/guess!
        lambda_r2=800u"nm", # TODO: This is a test/guess!
    )

    LNB_C = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full;
        miller_delta,
    )
    return LNB_C
end

export LNB_C
LNB_C = create_lnb_c()