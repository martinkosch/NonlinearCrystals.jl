function create_ktp_h()
    metadata = Dict(
        :description => "Hydrothermally-grown KTP (Potassium titanyl phosphate)",
        :formula => "KTiOP0₄",
        :pointgroup => "mm2",
        :lattice_params => (6.404u"Å", 10.616u"Å", 12.814u"Å"),
        :density => 2.945u"g/cm^3",
        :mohs_hardness => 5,
        :axes_assignment => (:X => :a, :Y => :b, :Z => :c),
    )

    # From: Handbook of Nonlinear Crystals
    lambda_min_bbo = 0.35u"µm"
    lambda_max_bbo = 4.5u"µm"

    n_x_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.1146 + 0.89188 * λ^2 / (λ^2 - (0.20861u"µm")^2) - 0.01320u"µm^-2" * λ^2),
        (lambda_min_bbo, lambda_max_bbo);
        dn_dT_fun=(λ, T) -> (0.1323u"µm^3" * λ^-3 - 0.4385u"µm^2" * λ^-2 + 1.2307u"µm" * λ^-1 + 0.7709) * 1e-5u"K^-1",
        T_ref=293.15u"K",
    )

    n_y_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.1518 + 0.87862 * λ^2 / (λ^2 - (0.21801u"µm")^2) - 0.01327u"µm^-2" * λ^2),
        (lambda_min_bbo, lambda_max_bbo);
        dn_dT_fun=(λ, T) -> (0.5014u"µm^3" * λ^-3 - 2.0030u"µm^2" * λ^-2 + 3.3016u"µm" * λ^-1 + 0.7498) * 1e-5u"K^-1",
        T_ref=293.15u"K",
    )

    n_z_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.3136 + 1.00012 * λ^2 / (λ^2 - (0.23831u"µm")^2) - 0.01679u"µm^-2" * λ^2),
        (lambda_min_bbo, lambda_max_bbo);
        dn_dT_fun=(λ, T) -> (0.3896u"µm^3" * λ^-3 - 1.3332u"µm^2" * λ^-2 + 2.2762u"µm" * λ^-1 + 2.1151) * 1e-5u"K^-1",
        T_ref=293.15u"K",
    )

    d = construct_d_tensor(metadata[:pointgroup]; d15=1.95u"pm/V", d24=3.9u"pm/V", d31=1.95u"pm/V", d32=3.9u"pm/V", d33=15.3u"pm/V")

    KTP_H = BidirectionalCrystal(
        metadata,
        n_x_principal,
        n_y_principal,
        n_z_principal,
        d,
    )
    return KTP_H
end

export KTP_H
KTP_H = create_ktp_h()