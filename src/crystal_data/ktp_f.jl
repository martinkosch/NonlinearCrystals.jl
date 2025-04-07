function create_ktp_f()
    metadata = Dict(
        :description => "Flux-grown KTP (Potassium titanyl phosphate)",
        :formula => "KTiOP0₄",
        :point_group => "mm2",
        :lattice_params => (6.404u"Å", 10.616u"Å", 12.814u"Å"),
        :density => 2.945u"g/cm^3",
        :mohs_hardness => 5,
    )

    # From: Handbook of Nonlinear Crystals
    lambda_min = 0.35u"µm"
    lambda_max = 4.5u"µm"
    temp_ref = 293.15u"K"

    n_X_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.0065 + 0.03901u"µm^2" / (λ^2 - 0.04251u"µm^2") - 0.01327u"µm^-2" * λ^2)
            dn_dtemp = (0.1323u"µm^3" * λ^-3 - 0.4385u"µm^2" * λ^-2 + 1.2307u"µm" * λ^-1 + 0.7709) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Y_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.0333 + 0.04154u"µm^2" / (λ^2 - 0.04547u"µm^2") - 0.01408u"µm^-2" * λ^2)
            dn_dtemp = (0.5014u"µm^3" * λ^-3 - 2.0030u"µm^2" * λ^-2 + 3.3016u"µm" * λ^-1 + 0.7498) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Z_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.3134 + 0.05694u"µm^2" / (λ^2 - 0.05658u"µm^2") - 0.01682u"µm^-2" * λ^2)
            dn_dtemp = (0.3896u"µm^3" * λ^-3 - 1.3332u"µm^2" * λ^-2 + 2.2762u"µm" * λ^-1 + 2.1151) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d31=1.95u"pm/V", d32=3.9u"pm/V", d33=15.3u"pm/V")
    # d_XYZ_full = calc_d_XYZ_full(metadata[:point_group], I(3), false; d31=3.7u"pm/V", d32=2.2u"pm/V", d33=14.6u"pm/V", d15=3.7u"pm/V", d24=1.9u"pm/V")
    
    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_X_principal, 
        n_Y_principal, 
        n_Z_principal, 
        temp_ref;  # TODO: This is a test/guess!
        lambda_r1=1064u"nm", # TODO: This is a test/guess!
        lambda_r2=1064u"nm", # TODO: This is a test/guess!
    )

    KTP_F = BidirectionalCrystal(
        metadata,
        n_X_principal,
        n_Y_principal,
        n_Z_principal,
        d_XYZ_full;
        miller_delta,
    )
    return KTP_F
end

export KTP_F
KTP_F = create_ktp_f()