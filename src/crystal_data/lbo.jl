function create_lbo()
    metadata = Dict(
        :description => "LBO (Lithium Triborate)",
        :formula => "LiB₃0₅",
        :point_group => "mm2",
        :lattice_params => (8.447u"Å", 7.380u"Å", 5.120u"Å"),
        :density => 2.47u"g/cm^3",
        :mohs_hardness => 6,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#kato94temperature
    lambda_min = 0.155u"µm"
    lambda_max = 3.2u"µm"
    temp_ref = 293.15u"K"

    n_X_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.4542 + 0.01125u"µm^2" / (λ^2 - 0.01135u"µm^2") - 0.01388u"µm^-2" * λ^2)
            delta_temp = T - temp_ref
            delta_n_temp = (-3.76u"µm^-1" * λ + 2.30) * 1e-6u"K^-1" * (delta_temp + 29.13e-3u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Y_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.5390 + 0.01277u"µm^2" / (λ^2 - 0.01189u"µm^2") - 0.01849u"µm^-2" * λ^2 + 4.3025e-5u"µm^-4" * λ^4 - 2.9131e-5u"µm^-6" * λ^6)
            delta_temp = T - temp_ref
            delta_n_temp = (6.01u"µm^-1" * λ - 19.40) * 1e-6u"K^-1" * (delta_temp - 32.89e-4u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    n_Z_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(2.5865 + 0.01310u"µm^2" / (λ^2 - 0.01223u"µm^2") - 0.01862u"µm^-2" * λ^2 + 4.5778e-5u"µm^-4" * λ^4 - 3.2526e-5u"µm^-6" * λ^6)
            delta_temp = T - temp_ref
            delta_n_temp = (1.50u"µm^-1" * λ - 9.70) * 1e-6u"K^-1" * (delta_temp - 74.49e-4u"K^-1" * delta_temp^2)
            return n_lam + delta_n_temp
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )

    rot_mat = rot_mat_crys_to_diel((:X, :Z, :Y))
    d_XYZ_full = calc_d_XYZ_full(metadata[:point_group], rot_mat; d31=-0.67u"pm/V", d32=0.85u"pm/V", d33=0.04u"pm/V")
    
    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_X_principal, 
        n_Y_principal, 
        n_Z_principal, 
        temp_ref;  # TODO: This is a test/guess!
        lambda_r1=800u"nm", # TODO: This is a test/guess!
        lambda_r2=800u"nm", # TODO: This is a test/guess!
    )

    LBO = BidirectionalCrystal(
        metadata,
        n_X_principal,
        n_Y_principal,
        n_Z_principal,
        d_XYZ_full;
        miller_delta,
    )
    return LBO
end

export LBO
const LBO = create_lbo()