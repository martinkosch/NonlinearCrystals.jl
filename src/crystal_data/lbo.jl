function create_lbo()
    metadata = Dict(
        :description => "LBO (Lithium Triborate)",
        :formula => "LiB₃0₅",
        :pointgroup => "mm2",
        :lattice_params => (8.447u"Å", 7.380u"Å", 5.120u"Å"),
        :density => 2.47u"g/cm^3",
        :mohs_hardness => 6,
        :axes_assignment => (:X => :a, :Y => :c, :Z => :b),
    )

    # From: Datasheet EKSMA OPTICS, temperature dependence from Handbook of Nonlinear Crystals
    lambda_min = 0.155u"µm"
    lambda_max = 3.2u"µm"

    n_x_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.4542 + 0.01125u"µm^2" / (λ^2 - 0.01135u"µm^2") - 0.01388u"µm^-2" * λ^2),
        (lambda_min, lambda_max);
        dn_dtemp_fun=(λ,T) -> -9.3e-6u"K^-1",
        temp_ref=293.15u"K",
    )

    n_y_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.5390 + 0.01277u"µm^2" / (λ^2 - 0.01189u"µm^2") - 0.01849u"µm^-2" * λ^2 + 4.3025e-5u"µm^-4" * λ^4 - 2.9131e-5u"µm^-6" * λ^6),
        (lambda_min, lambda_max);
        dn_dtemp_fun=(λ,T) -> -13.6e-6u"K^-1",
        temp_ref=293.15u"K",
    )

    n_z_principal = SellmeierFunction(
        (λ, T) -> sqrt(2.5865 + 0.01310u"µm^2" / (λ^2 - 0.01223u"µm^2") - 0.01862u"µm^-2" * λ^2 + 4.5778e-5u"µm^-4" * λ^4 - 3.2526e-5u"µm^-6" * λ^6),
        (lambda_min, lambda_max);
        dn_dtemp_fun=(λ,T) -> (-6.3 - 2.1u"µm^-1" * λ) * 1e-6u"K^-1",
        temp_ref=293.15u"K",
    )

    d = construct_d_tensor(metadata[:pointgroup]; d31=1.05u"pm/V", d32=-0.98u"pm/V", d33=0.05u"pm/V") 

    LBO = BidirectionalCrystal(
        metadata,
        n_x_principal,
        n_y_principal,
        n_z_principal,
        d,
    )
    return LBO
end

export LBO
LBO = create_lbo()