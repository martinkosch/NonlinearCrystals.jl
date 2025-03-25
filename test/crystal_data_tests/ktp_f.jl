using NonlinearCrystals
using Test
using Unitful

@test KTP_F.n_x_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1
@test KTP_F.n_y_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1
@test KTP_F.n_z_principal(default_lambda(KTP_F), default_temp(KTP_F)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(KTP_F.n_x_principal(0.4047u"µm", 293u"K"), 1.8249, atol=0.01)
@test isapprox(KTP_F.n_x_principal(0.5893u"µm", 293u"K"), 1.7689, atol=0.01)
@test isapprox(KTP_F.n_x_principal(1.3414u"µm", 293u"K"), 1.7314, atol=0.01)

@test isapprox(KTP_F.n_y_principal(0.4047u"µm", 293u"K"), 1.8410, atol=0.01)
@test isapprox(KTP_F.n_y_principal(0.5893u"µm", 293u"K"), 1.7780, atol=0.01)
@test isapprox(KTP_F.n_y_principal(1.3414u"µm", 293u"K"), 1.7387, atol=0.01)

@test isapprox(KTP_F.n_z_principal(0.4047u"µm", 293u"K"), 1.9629, atol=0.01)
@test isapprox(KTP_F.n_z_principal(0.5893u"µm", 293u"K"), 1.8740, atol=0.01)
@test isapprox(KTP_F.n_z_principal(1.3414u"µm", 293u"K"), 1.8211, atol=0.01)

# Test optical axis
@test isapprox(optical_axis_angle(KTP_F, 0.5461u"µm"), 37.4u"°" / 2, atol=ustrip(u"rad", 1u"°"))

# Test temperature derivatives
@test isapprox(KTP_F.n_x_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 2.41e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_F.n_x_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 1.65e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(KTP_F.n_y_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 3.21e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_F.n_y_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 2.50e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(KTP_F.n_z_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 4.27e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_F.n_z_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 3.40e-5u"K^-1", atol=0.1e-5u"K^-1")

# Test sampled phasematches (from: Handbook of Nonlinear Crystals)
pm1 = find_nearest_pm_along_theta_phi(90u"°", 24.59u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=1.0642u"µm", lambda_b=0.5321u"µm", temp=293u"K") 
@test all(isnothing.(pm1.o_or_e_r1_r2_b)) || all(pm1.o_or_e_r1_r2_b .== [:e, :o, :e])
@test isapprox(pm1.phi_pm, 24.59u"°", atol=ustrip(u"rad", 2u"°"))
@test isapprox(pm1.walkoff_angle_r1_r2_b[1], 0.202u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm1.walkoff_angle_r1_r2_b[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm1.walkoff_angle_r1_r2_b[3], 0.268u"°", atol=ustrip(u"rad", 0.1u"°")) 

pm2 = find_nearest_pm_along_theta_phi(57.95u"°", 90.0u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=2.9365u"µm", lambda_b=2.9365u"µm" / 2, temp=293u"K")
@test all(isnothing.(pm2.o_or_e_r1_r2_b)) || all(pm2.o_or_e_r1_r2_b .== [:o, :e, :o])
@test isapprox(pm2.theta_pm, 57.95u"°", atol=ustrip(u"rad", 4u"°"))

pm3 = find_nearest_pm_along_theta_phi(56.22u"°", 90.0u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=1.2u"µm", lambda_b=1.2u"µm" / 2, temp=293u"K")
@test all(isnothing.(pm3.o_or_e_r1_r2_b)) || all(pm3.o_or_e_r1_r2_b .== [:o, :e, :o])
@test isapprox(pm3.theta_pm, 56.22u"°", atol=ustrip(u"rad", 4u"°"))

pm4 = find_nearest_pm_along_theta_phi(67.47u"°", 0.0u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=1.2u"µm", lambda_b=1.2u"µm" / 2, temp=293u"K")
@test all(isnothing.(pm4.o_or_e_r1_r2_b)) || all(pm4.o_or_e_r1_r2_b .== [:o, :e, :o])
@test isapprox(pm4.theta_pm, 67.47u"°", atol=ustrip(u"rad", 4u"°"))
