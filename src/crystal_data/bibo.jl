function create_bibo()
    metadata = Dict(
        :description => "BIBO (Bismuth triborate)",
        :formula => "BiB₃O₆",
        :point_group => "2",
        :lattice_params => (7.116u"Å", 4.993u"Å", 6.508u"Å"),
        :density => 5.033u"g/cm^3",
        :mohs_hardness => 5,
    )

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#miyata09
    temp_ref = 293.15u"K"
    lambda_min = 0.3263u"µm"
    lambda_max = 3.083u"µm"
    
    n_X_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.0759 + 0.03169u"µm^2" / (λ^2 - 0.03323u"µm^2") - 0.01402u"µm^-2" * λ^2)
            dndT = (0.1178u"µm^3" * λ^-3 - 0.2282u"µm^2" * λ^-2 + 0.7965u"µm" * λ^-1 - 0.2962) * 1e-5u"K^-1"
            return n_lam + dndT * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )
    
    n_Y_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.1698 + 0.03666u"µm^2" / (λ^2 - 0.03599u"µm^2") - 0.01819u"µm^-2" * λ^2)
            dndT = (0.1777u"µm^3" * λ^-3 - 0.5594u"µm^2" * λ^-2 + 0.8336u"µm" * λ^-1 - 0.9960) * 1e-5u"K^-1"
            return n_lam + dndT * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )
    
    n_Z_principal = SellmeierFunction(
        (λ, T) -> begin
            n_lam = sqrt(3.6546 + 0.05116u"µm^2" / (λ^2 - 0.03713u"µm^2") - 0.02299u"µm^-2" * λ^2)
            dndT = (0.2075u"µm^3" * λ^-3 - 0.7026u"µm^2" * λ^-2 + 0.8378u"µm" * λ^-1 - 1.0210) * 1e-5u"K^-1"
            return n_lam + dndT * (T - temp_ref)
        end,
        (lambda_min, lambda_max);
        temp_ref,
    )    

    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#hellwig00linear
    phi = 47u"°"
    rot_mat = rot_mat_crys_to_diel((:Z, :X, :Y); rotate_about=:X, phi)
    d_XYZ_full = calc_d_XYZ_full(
        metadata[:point_group],
        rot_mat,
        false;
        d21=2.32u"pm/V", # (2, 1, 1)
        d16=2.8u"pm/V", # (1, 1, 2)
        d22=2.538u"pm/V", # (2, 2, 2)
        d13=0.0u"pm/V", # (1, 3, 3) 
        d23=-1.31u"pm/V", # (2, 3, 3)
        d34=-0.91u"pm/V", # (3, 3, 2)
        d14=2.43u"pm/V", # (1, 3, 2)
        d25=2.32u"pm/V", # (2, 3, 1)
        d36=2.43u"pm/V", # (3, 1, 2)
    )

    miller_delta = calc_miller_delta(
        d_XYZ_full, 
        n_X_principal, 
        n_Y_principal, 
        n_Z_principal, 
        temp_ref;  
        lambda_r1=1079.5u"nm", 
        lambda_r2=1079.5u"nm", 
    )

    BIBO = BidirectionalCrystal(
        metadata,
        n_X_principal,
        n_Y_principal,
        n_Z_principal,
        d_XYZ_full;
        miller_delta,
    )
    return BIBO
end

export BIBO
BIBO = create_bibo()