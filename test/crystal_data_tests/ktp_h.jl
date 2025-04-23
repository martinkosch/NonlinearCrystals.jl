using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_ktp_h()) <: BidirectionalCrystal
@test isa(crystal_system(KTP_H), String)

@test KTP_H.n_X_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1
@test KTP_H.n_Y_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1
@test KTP_H.n_Z_principal(default_lambda(KTP_H), default_temp(KTP_H)) > 1

# Sample refractive indices 
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
@test isapprox(KTP_H.n_X_principal(0.53u"µm", 293u"K"), 1.7787, atol=0.01) 
@test isapprox(KTP_H.n_X_principal(1.06u"µm", 293u"K"), 1.7400, atol=0.01) 

@test isapprox(KTP_H.n_Y_principal(0.53u"µm", 293u"K"), 1.7924, atol=0.01) 
@test isapprox(KTP_H.n_Y_principal(1.06u"µm", 293u"K"), 1.7469, atol=0.01) 

@test isapprox(KTP_H.n_Z_principal(0.53u"µm", 293u"K"), 1.8873, atol=0.01) 
@test isapprox(KTP_H.n_Z_principal(1.06u"µm", 293u"K"), 1.8304, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(KTP_H, 0.5461u"µm"), 37.4u"°"/2, atol=ustrip(u"rad", 1u"°")) 

# Test temperature derivatives
@test isapprox(dn_dtemp(KTP_H.n_X_principal, 0.5321u"µm", 293u"K"), 2.41e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_H.n_X_principal, 1.0642u"µm", 293u"K"), 1.65e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(dn_dtemp(KTP_H.n_Y_principal, 0.5321u"µm", 293u"K"), 3.21e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_H.n_Y_principal, 1.0642u"µm", 293u"K"), 2.50e-5u"K^-1", atol=0.1e-5u"K^-1")

@test isapprox(dn_dtemp(KTP_H.n_Z_principal, 0.5321u"µm", 293u"K"), 4.27e-5u"K^-1", atol=0.1e-5u"K^-1")
@test isapprox(dn_dtemp(KTP_H.n_Z_principal, 1.0642u"µm", 293u"K"), 3.40e-5u"K^-1", atol=0.1e-5u"K^-1")

# Test sampled phasematches
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
pm1 = find_nearest_pm_along_theta_phi(90u"°", 24.0u"°", (:lo, :hi, :lo), KTP_H; lambda_r1=1.062u"µm", lambda_b=1.062u"µm" / 2, temp=293u"K") 
@test all(pm1.pm_type[1].o_or_e_rrb .== (:e, :o, :e))
@test pm1.pm_type[1].principal_plane == :XY
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.phi_pm, 24.0u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm1.eff_data.d_eff), 3.48u"pm/V", rtol=0.2) # From SNLO

pm2 = find_nearest_pm_along_theta_phi(63.2u"°", 90.0u"°", (:lo, :hi, :lo), KTP_H; lambda_r1=1.338u"µm", lambda_b=0.446u"µm", temp=293u"K")
@test all(pm2.pm_type[1].o_or_e_rrb .== (:o, :e, :o))
@test pm2.pm_type[1].principal_plane == :YZ
@test isnothing(pm2.pm_type[2])
@test isapprox(pm2.theta_pm, 63.2u"°", atol=ustrip(u"rad", 4u"°"))
@test isapprox(abs(pm2.eff_data.d_eff), 1.81u"pm/V", rtol=0.2) # From SNLO

# Test sampled noncritical phasematches
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
pm3s = find_all_ncpm_over_lambda((:lo,:hi,:lo), KTP_H, 47u"°C"; lambda_b=0.4396u"µm", principal_axis=:X)
@test length(pm3s) == 1
pm3 = pm3s[1]
@test all(isapprox(abs.(pm3.E_dir_rrb[1]), [0, 1, 0]))
@test all(isapprox(abs.(pm3.E_dir_rrb[2]), [0, 0, 1]))
@test all(isapprox(abs.(pm3.E_dir_rrb[3]), [0, 1, 0]))
@test pm3.pm_type[1].type == pm3.pm_type[2].type == "II/III"
@test isapprox(pm3.theta_pm, 90.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm3.phi_pm, 0.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm3.lambda_rrb[1], 1318.8u"nm", atol=50u"nm")
@test isapprox(pm3.bw_data.temp_L_bw, 8.5u"K * cm", atol=5u"K * cm")

pm4s = find_all_ncpm_over_lambda((:lo,:hi,:lo), KTP_H, 463u"°C"; lambda_b=0.446u"µm", principal_axis=:X)
@test length(pm4s) == 1
pm4 = pm4s[1]
@test all(isapprox(abs.(pm4.E_dir_rrb[1]), [0, 1, 0]))
@test all(isapprox(abs.(pm4.E_dir_rrb[2]), [0, 0, 1]))
@test all(isapprox(abs.(pm4.E_dir_rrb[3]), [0, 1, 0]))
@test pm4.pm_type[1].type == pm4.pm_type[2].type == "II/III"
@test isapprox(pm4.theta_pm, 90.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm4.phi_pm, 0.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm4.lambda_rrb[1], 1338.0u"nm", atol=50u"nm")
@test isapprox(pm4.bw_data.temp_L_bw, 8.5u"K * cm", atol=5u"K * cm")

pm5s = find_all_ncpm_over_lambda((:lo,:hi,:lo), KTP_H, 20u"°C"; lambda_b=0.49715u"µm", principal_axis=:Y)
@test length(pm5s) == 1
pm5 = pm5s[1]
@test all(isapprox(abs.(pm5.E_dir_rrb[1]), [1, 0, 0]))
@test all(isapprox(abs.(pm5.E_dir_rrb[2]), [0, 0, 1]))
@test all(isapprox(abs.(pm5.E_dir_rrb[3]), [1, 0, 0]))
@test pm5.pm_type[1].type == pm5.pm_type[2].type == "II/III"
@test isapprox(pm5.theta_pm, 90.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm5.phi_pm, 90.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm5.lambda_rrb[1], 994.3u"nm", atol=50u"nm")
@test isapprox(pm5.bw_data.temp_L_bw, 175u"K * cm", atol=5u"K * cm")
