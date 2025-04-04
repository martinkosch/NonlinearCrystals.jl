using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_ktp_f()) <: BidirectionalCrystal
@test isa(crystal_system(KTP_F), String)

@test KTP_F.n_X_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1
@test KTP_F.n_Y_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1
@test KTP_F.n_Z_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(KTP_F.n_X_principal(0.4047u"µm", 293u"K"), 1.8249, atol=0.01)
@test isapprox(KTP_F.n_X_principal(0.5893u"µm", 293u"K"), 1.7689, atol=0.01)
@test isapprox(KTP_F.n_X_principal(1.3414u"µm", 293u"K"), 1.7314, atol=0.01)

@test isapprox(KTP_F.n_Y_principal(0.4047u"µm", 293u"K"), 1.8410, atol=0.01)
@test isapprox(KTP_F.n_Y_principal(0.5893u"µm", 293u"K"), 1.7780, atol=0.01)
@test isapprox(KTP_F.n_Y_principal(1.3414u"µm", 293u"K"), 1.7387, atol=0.01)

@test isapprox(KTP_F.n_Z_principal(0.4047u"µm", 293u"K"), 1.9629, atol=0.01)
@test isapprox(KTP_F.n_Z_principal(0.5893u"µm", 293u"K"), 1.8740, atol=0.01)
@test isapprox(KTP_F.n_Z_principal(1.3414u"µm", 293u"K"), 1.8211, atol=0.01)

# Test optical axis
@test isapprox(optical_axis_angle(KTP_F, 0.5461u"µm"), 37.4u"°" / 2, atol=ustrip(u"rad", 1u"°"))

# Test temperature derivatives
@test isapprox(dn_dtemp(KTP_F.n_X_principal, 0.5321u"µm", 293u"K"), 2.41e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_F.n_X_principal, 1.0642u"µm", 293u"K"), 1.65e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(dn_dtemp(KTP_F.n_Y_principal, 0.5321u"µm", 293u"K"), 3.21e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_F.n_Y_principal, 1.0642u"µm", 293u"K"), 2.50e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(dn_dtemp(KTP_F.n_Z_principal, 0.5321u"µm", 293u"K"), 4.27e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_F.n_Z_principal, 1.0642u"µm", 293u"K"), 3.40e-5u"K^-1", atol=0.1e-5u"K^-1")

# Test sampled phasematches (from: Handbook of Nonlinear Crystals)
pm1 = find_nearest_pm_along_theta_phi(90u"°", 24.59u"°", (:lo, :hi, :lo), KTP_F; lambda_r1=1.0642u"µm", lambda_b=0.5321u"µm", temp=293u"K") 
@test all(pm1.pm_type[1].o_or_e_rrb .== (:e, :o, :e))
@test pm1.pm_type[1].principal_plane == :XY
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.phi_pm, 24.59u"°", atol=ustrip(u"rad", 2u"°"))
@test isapprox(pm1.walkoff_angle_rrb[1], 0.202u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[3], 0.268u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(abs(pm1.eff_data.d_eff), 3.59u"pm/V", rtol=0.2) # From SNLO

pm2 = find_nearest_pm_along_theta_phi(34.5u"°", 90.0u"°", (:hi, :hi, :lo), KTP_F; lambda_r1=2.9365u"µm", lambda_b=2.9365u"µm" / 2, temp=293u"K")
@test all(pm2.pm_type[1].o_or_e_rrb .== (:e, :e, :o))
@test pm2.pm_type[1].principal_plane == :YZ
@test isnothing(pm2.pm_type[2])
@test isapprox(pm2.theta_pm, 34.5u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm2.eff_data.d_eff), 0.0u"pm/V", rtol=0.2) # From SNLO

pm3 = find_nearest_pm_along_theta_phi(57.95u"°", 90.0u"°", (:lo, :hi, :lo), KTP_F; lambda_r1=2.9365u"µm", lambda_b=2.9365u"µm" / 2, temp=293u"K")
@test all(pm3.pm_type[1].o_or_e_rrb .== (:o, :e, :o))
@test pm3.pm_type[1].principal_plane == :YZ
@test isnothing(pm3.pm_type[2])
@test isapprox(pm3.theta_pm, 57.95u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm3.eff_data.d_eff), 1.29u"pm/V", rtol=0.3) # From SNLO

pm4 = find_nearest_pm_along_theta_phi(56.22u"°", 90.0u"°", (:lo, :hi, :lo), KTP_F; lambda_r1=1.2u"µm", lambda_b=1.2u"µm" / 2, temp=293u"K")
@test all(pm4.pm_type[1].o_or_e_rrb .== (:o, :e, :o))
@test pm4.pm_type[1].principal_plane == :YZ
@test isnothing(pm4.pm_type[2])
@test isapprox(pm4.theta_pm, 56.22u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm4.eff_data.d_eff), 1.54u"pm/V", rtol=0.2) # From SNLO

pm5 = find_nearest_pm_along_theta_phi(67.47u"°", 0.0u"°", (:lo, :hi, :lo), KTP_F; lambda_r1=1.2u"µm", lambda_b=1.2u"µm" / 2, temp=293u"K")
@test all(pm5.pm_type[1].o_or_e_rrb .== (:o, :e, :o))
@test pm5.pm_type[1].principal_plane == :XZ
@test isnothing(pm5.pm_type[2])
@test isapprox(pm5.theta_pm, 67.47u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm5.eff_data.d_eff), 3.42u"pm/V", rtol=0.2) # From SNLO