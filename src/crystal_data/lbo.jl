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

    # From: Datasheet EKSMA OPTICS, temperature dependence from Kato, “Temperature-tuned 90º phase-matching properties of LiB3O5,” IEEE J. 1994
    lambda_min = 0.155u"µm"
    lambda_max = 3.2u"µm"
    temp_ref = 293.15u"K"

    n_x_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.4542 + 0.01125u"µm^2" / (λ^2 - 0.01135u"µm^2") - 0.01388u"µm^-2" * λ^2)
            delta_temp = T - temp_ref
            delta_n_temp = (-3.76u"µm^-1" * λ + 2.30) * 1e-6u"K^-1" * (delta_temp + 29.13e-3u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_y_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.5390 + 0.01277u"µm^2" / (λ^2 - 0.01189u"µm^2") - 0.01849u"µm^-2" * λ^2 + 4.3025e-5u"µm^-4" * λ^4 - 2.9131e-5u"µm^-6" * λ^6)
            delta_temp = T - temp_ref
            delta_n_temp = (6.01u"µm^-1" * λ - 19.40) * 1e-6u"K^-1" * (delta_temp - 32.89e-4u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_z_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.5865 + 0.01310u"µm^2" / (λ^2 - 0.01223u"µm^2") - 0.01862u"µm^-2" * λ^2 + 4.5778e-5u"µm^-4" * λ^4 - 3.2526e-5u"µm^-6" * λ^6)
            delta_temp = T - temp_ref
            delta_n_temp = (1.50u"µm^-1" * λ - 9.70) * 1e-6u"K^-1" * (delta_temp - 74.49e-4u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
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