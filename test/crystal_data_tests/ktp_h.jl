using NonlinearCrystals
using Test
using Unitful

@test KTP_H.n_x_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1
@test KTP_H.n_y_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1
@test KTP_H.n_z_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(KTP_H.n_x_principal(0.53u"µm", 293u"K"), 1.7787, atol=0.01) 
@test isapprox(KTP_H.n_x_principal(1.06u"µm", 293u"K"), 1.7400, atol=0.01) 

@test isapprox(KTP_H.n_y_principal(0.53u"µm", 293u"K"), 1.7924, atol=0.01) 
@test isapprox(KTP_H.n_y_principal(1.06u"µm", 293u"K"), 1.7469, atol=0.01) 

@test isapprox(KTP_H.n_z_principal(0.53u"µm", 293u"K"), 1.8873, atol=0.01) 
@test isapprox(KTP_H.n_z_principal(1.06u"µm", 293u"K"), 1.8304, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(KTP_H, 0.5461u"µm"), 37.4u"°"/2, atol=ustrip(u"rad", 1u"°")) 

# Test temperature derivatives
@test isapprox(KTP_H.n_x_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 2.41e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_H.n_x_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 1.65e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(KTP_H.n_y_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 3.21e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_H.n_y_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 2.50e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(KTP_H.n_z_principal.dn_dtemp_fun(0.5321u"µm", 293u"K"), 4.27e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(KTP_H.n_z_principal.dn_dtemp_fun(1.0642u"µm", 293u"K"), 3.40e-5u"K^-1", atol=0.1e-5u"K^-1")

# Test sampled phasematches (from: Handbook of Nonlinear Crystals)
pm1 = find_nearest_pm_along_theta_phi(90u"°", 24.0u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=1.062u"µm", lambda_b=1.062u"µm" / 2, temp=293u"K") 
@test all(isnothing.(pm1.o_or_e_r1_r2_b)) || all(pm1.o_or_e_r1_r2_b .== [:e, :o, :e])
@test isapprox(pm1.phi_pm, 24.0u"°", atol=ustrip(u"rad", 4u"°"))

pm2 = find_nearest_pm_along_theta_phi(63.2u"°", 90.0u"°", [:lo, :hi, :lo], KTP_H; lambda_r1=1.338u"µm", lambda_b=0.446u"µm", temp=293u"K")
@test all(isnothing.(pm2.o_or_e_r1_r2_b)) || all(pm2.o_or_e_r1_r2_b .== [:o, :e, :o])
@test isapprox(pm2.theta_pm, 63.2u"°", atol=ustrip(u"rad", 4u"°"))
