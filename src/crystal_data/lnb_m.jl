# TODO: Search for proper source with complete values
# function create_lnb_m()
#     metadata = Dict(
#         :description => "LNB (Lithium Niobate, Magnesium-Oxide-Doped with 5 mole % MgO)",
#         :formula => "MgO:LiNbO₃",
#         :point_group => "3m",
#         :lattice_params => (5.148u"Å", 5.148u"Å", 3.863u"Å"),
#         :density => 4.628u"g/cm^3",
#         :mohs_hardness => 5,
#     )

#     # Handbook of Nonlinear Crystals
#     lambda_min = 0.4u"µm"
#     lambda_max = 5.5u"µm"

#     n_o_principal = SellmeierFunction(
#         (λ, T) -> sqrt(
#             4.9017 +
#             0.112280 / (λ^2 - 0.049656u"µm^2") -
#             0.039636u"µm^-2" * λ^2
#         ),
#         (lambda_min, lambda_max);
#         temp_ref = 293.15u"K",
#     )

#     n_e_principal = SellmeierFunction(
#         (λ, T) -> sqrt(
#             4.5583 +
#             0.091806 / (λ^2 - 0.048086u"µm^2") -
#             0.032068u"µm^-2" * λ^2
#         ),
#         (lambda_min, lambda_max);
#         temp_ref = 293.15u"K",
#     )
    

#     d_XYZ_full = calc_d_XYZ_full(metadata[:point_group]; d22=2.46u"pm/V", d31=-4.64u"pm/V", d33=-41.7u"pm/V") # Measured at 1.058 µm

# miller_delta = calc_miller_delta(
#     d_XYZ_full, 
#     n_o_principal,
#     n_e_principal,
#     temp_ref;  # TODO: This is a test/guess!
#     lambda_r1=800u"nm", # TODO: This is a test/guess!
#     lambda_r2=800u"nm", # TODO: This is a test/guess!
# )

#     LNB_M = UnidirectionalCrystal(
#         metadata,
#         n_o_principal,
#         n_e_principal,
#         d_XYZ_full;
        # miller_delta,
#     )
#     return LNB_M
# end

# export LNB_M
# LNB_M = create_lnb_m()