function create_ktp_h()
    metadata = Dict(
        :description => "Hydrothermally-grown KTP (Potassium titanyl phosphate)",
        :formula => "KTiOP0₄",
        :pointgroup => "mm2",
        :lattice_params => (6.404u"Å", 10.616u"Å", 12.814u"Å"),
        :density => 2.945u"g/cm^3",
        :mohs_hardness => 5,
        :axes_assignment_XYZ => (:a, :b, :c)
    )

    # From: Handbook of Nonlinear Crystals
    lambda_min = 0.35u"µm"
    lambda_max = 4.5u"µm"
    temp_ref = 293.15u"K"

    n_X_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.1146 + 0.89188 * λ^2 / (λ^2 - (0.20861u"µm")^2) - 0.01320u"µm^-2" * λ^2)
            dn_dtemp = (0.1323u"µm^3" * λ^-3 - 0.4385u"µm^2" * λ^-2 + 1.2307u"µm" * λ^-1 + 0.7709) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Y_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.1518 + 0.87862 * λ^2 / (λ^2 - (0.21801u"µm")^2) - 0.01327u"µm^-2" * λ^2)
            dn_dtemp = (0.5014u"µm^3" * λ^-3 - 2.0030u"µm^2" * λ^-2 + 3.3016u"µm" * λ^-1 + 0.7498) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Z_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.3136 + 1.00012 * λ^2 / (λ^2 - (0.23831u"µm")^2) - 0.01679u"µm^-2" * λ^2)
            dn_dtemp = (0.3896u"µm^3" * λ^-3 - 1.3332u"µm^2" * λ^-2 + 2.2762u"µm" * λ^-1 + 2.1151) * 1e-5u"K^-1"
            return n_lam + dn_dtemp * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    d_XYZ_full = calc_d_XYZ_full(metadata[:pointgroup]; d31=1.95u"pm/V", d32=3.9u"pm/V", d33=15.3u"pm/V")


    KTP_H = BidirectionalCrystal(
        metadata,
        n_X_principal,
        n_Y_principal,
        n_Z_principal,
        d_XYZ_full,
    )
    return KTP_H
end

export KTP_H
KTP_H = create_ktp_h()