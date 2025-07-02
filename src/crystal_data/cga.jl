function create_cga()
    metadata = Dict(
        :description => "CGA (Cadmium germanium arsenide)",
        :formula => "CdGeAs₂",
        :point_group => "-42m",
        :lattice_params => (1.26512u"nm", 1.26512u"nm", 2.54435u"nm"),
        :density => 5.6u"g/cm^3",
        :mohs_hardness => 4,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    lambda_min = 2.4u"µm"
    lambda_max = 18.0u"µm"
    temp_ref = 293u"K"

    n_o_principal = SellmeierFunction(
        (λ, T) -> sqrt(1 + 9.1064 + 2.2998u"µm^2" / (λ^2 - 1.0872u"µm^2") + 1.62478u"µm^2" / (λ^2 - 1370u"µm^2")),
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_e_principal = SellmeierFunction(
        (λ, T) -> sqrt(1 + 10.8018 + 1.2152u"µm^2" / (λ^2 - 2.6971u"µm^2") + 1.6922u"µm^2" / (λ^2 - 1370u"µm^2")),
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d14=235.0u"pm/V")

    CGA = UnidirectionalCrystal(
        metadata,
        n_o_principal,
        n_e_principal,
        d_XYZ_full,
    )
    return CGA
end

export CGA
const CGA = create_cga()