function create_lnb_s()
    metadata = Dict(
        :description => "LNB (Lithium Niobate, stoichiometric melt with mole ratio Li/Nb = 1.000)",
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
            4.9130 +
            (0.1173u"µm^2" + 1.65e-8u"µm^2 * K^-2" * T^2) /
            (λ^2 - (0.212u"µm" + 2.7e-8u"µm * K^-2" * T^2)^2) -
            0.0278u"µm^-2" * λ^2
        ),
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_e_principal = SellmeierFunction(
        (λ, T) -> sqrt(
            4.5567 +
            2.605e-7u"K^-2" * T^2 +
            (0.0970u"µm^2" + 2.70e-8u"µm^2 * K^-2" * T^2) /
            (λ^2 - (0.201u"µm" + 5.4e-8u"µm * K^-2" * T^2)^2) -
            0.0224u"µm^-2" * λ^2
        ),
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d22=2.46u"pm/V", d31=-4.64u"pm/V", d33=-41.7u"pm/V") # Measured at 1.058 µm

    miller_delta = calc_miller_delta(
        d_XYZ_full,
        n_o_principal,
        n_e_principal,
        temp_ref;  # TODO: This is a test/guess!
        lambda_r1=800u"nm", # TODO: This is a test/guess!
        lambda_r2=800u"nm", # TODO: This is a test/guess!
    )

    LNB_S = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full;
        miller_delta,
    )
    return LNB_S
end

export LNB_S
LNB_S = create_lnb_s()