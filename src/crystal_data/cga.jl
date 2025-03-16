function create_cga()
    metadata = Dict(
        :description => "CGA (Cadmium germanium arsenide)",
        :formula => "CdGeAs₂",
        :pointgroup => "-42m",
        :lattice_params => (1.26512u"nm", 1.26512u"nm", 2.54435u"nm"),
        :density => 5.6u"g/cm^3",
        :mohs_hardness => 4,
    )

    # From Tamošauskas, Gintaras; Beresnevičius, Gvidas; Gadonas, Darius; Dubietis, Audrius . (2018). Transmittance and phase matching of BBO crystal in the 3−5 μm range and its application for the characterization of mid-infrared laser pulses. Optical Materials Express, 8(6), 1410–. doi:10.1364/ome.8.001410 
    lambda_min = 2.4u"µm"
    lambda_max = 18.0u"µm"

    n_o_principal = SellmeierFunction(
        # [9.1064, 2.2998, 1.6247],
        # [0.0, 1.0872, 1370] * 1u"µm^2",
        (λ,T) -> sqrt(1 + 9.1064 + 2.2998u"µm^2" / (λ^2 - 1.0872u"µm^2") + 1.62478u"µm^2" / (λ^2 - 1370u"µm^2")),
        (lambda_min, lambda_max);
        dn_dT_fun=(λ,T)->0.0u"K^-1",
        T_ref=293u"K",
    )

    n_e_principal = SellmeierFunction(
        # [10.8018, 1.2152, 1.6922],
        # [0.0, 2.6971, 1370] * 1u"µm^2",
        (λ,T) -> sqrt(1 + 10.8018 + 1.2152u"µm^2" / (λ^2 - 2.6971u"µm^2") + 1.6922u"µm^2" / (λ^2 - 1370u"µm^2")),
        (lambda_min, lambda_max);
        dn_dT_fun=(λ,T)->0.0u"K^-1",
        T_ref=293u"K",
    )

    d = [0 0 0 235.0 0 0;
        0 0 0 0 235.0 0;
        0 0 0 0 0 235.0] * 1u"pm/V"

    CGA = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d,
    )
    return CGA
end

export CGA
CGA = create_cga()