function create_bbo()
    metadata = Dict(
        :description => "BBO (Beta-barium borate)",
        :formula => "β-BaB₂O₄",
        :point_group => "3m",
        :lattice_params => (12.532u"Å", 12.532u"Å", 12.717u"Å"),
        :density => 3.85u"g/cm^3",
        :mohs_hardness => 4,
    )

    # From Tamošauskas, Gintaras; Beresnevičius, Gvidas; Gadonas, Darius; Dubietis, Audrius . (2018). Transmittance and phase matching of BBO crystal in the 3−5 μm range and its application for the characterization of mid-infrared laser pulses. Optical Materials Express, 8(6), 1410–. doi:10.1364/ome.8.001410 
    # lambda_min_bbo = 0.188u"µm"
    # lambda_max_bbo = 5.2u"µm"

    # n_o_bbo = SellmeierCoefficients(
    #     [0.90291, 0.83155, 0.76536], 
    #     [0.003926, 0.018786, 60.01] * 1u"µm^2", 
    #     (lambda_min_bbo, lambda_max_bbo); 
    #     dn_dT = -16.6e-6u"K^-1",
    # )

    # n_e_bbo = SellmeierCoefficients(
    #     [1.151075, 0.21803, 0.656], 
    #     [0.007142, 0.02259, 263] * 1u"µm^2", 
    #     (lambda_min_bbo, lambda_max_bbo); 
    #     dn_dT = -9.3e-6u"K^-1"
    # )

    # Handbook of Nonlinear Crystals
    lambda_min = 0.189u"µm"
    lambda_max = 3.5u"µm"
    temp_ref=293.15u"K"

    n_o_principal = SellmeierFunction(
        (λ, T) -> begin
            sqrt(2.7359 + 0.01878u"µm^2" / (λ^2 - 0.01822u"µm^2") - 0.01354u"µm^-2" * λ^2) - 16.6e-6u"K^-1" * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_e_principal = SellmeierFunction(
        (λ, T) -> begin
            sqrt(2.3753 + 0.01224u"µm^2" / (λ^2 - 0.01667u"µm^2") - 0.01516u"µm^-2" * λ^2) - 9.3e-6u"K^-1" * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d22=-2.2u"pm/V", d15=0.08u"pm/V", d33=0.0u"pm/V") 

    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_o_principal, 
        n_e_principal, 
        temp_ref; 
        lambda_r1=1064u"nm",
        lambda_r2=1064u"nm",
    )

    BBO = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full;
        miller_delta,
    )
    return BBO
end

export BBO
BBO = create_bbo()