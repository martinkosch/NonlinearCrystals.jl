function create_bbo()
    metadata = Dict(
        :description => "BBO (Beta-barium borate)",
        :formula => "β-BaB₂O₄",
        :point_group => "3m",
        :lattice_params => (12.532u"Å", 12.532u"Å", 12.717u"Å"),
        :density => 3.85u"g/cm^3",
        :mohs_hardness => 4,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    lambda_min = 0.189u"µm"
    lambda_max = 3.5u"µm"
    temp_ref = 293.15u"K"

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