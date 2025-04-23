using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_lbo()) <: BidirectionalCrystal
@test isa(crystal_system(LBO), String)

@test LBO.n_X_principal(default_lambda(LBO), default_temp(LBO)) > 1
@test LBO.n_Y_principal(default_lambda(LBO), default_temp(LBO)) > 1
@test LBO.n_Z_principal(default_lambda(LBO), default_temp(LBO)) > 1

# Sample refractive indices
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
@test isapprox(LBO.n_X_principal(0.2537u"µm", 293u"K"), 1.6335, atol=0.01) 
@test isapprox(LBO.n_X_principal(0.5321u"µm", 293u"K"), 1.57868, atol=0.01) 
@test isapprox(LBO.n_X_principal(1.1000u"µm", 293u"K"), 1.56432, atol=0.01) 

@test isapprox(LBO.n_Y_principal(0.2537u"µm", 293u"K"), 1.6582, atol=0.01) 
@test isapprox(LBO.n_Y_principal(0.5321u"µm", 293u"K"), 1.60642, atol=0.01) 
@test isapprox(LBO.n_Y_principal(1.1000u"µm", 293u"K"), 1.59005, atol=0.01) 

@test isapprox(LBO.n_Z_principal(0.2537u"µm", 293u"K"), 1.6792, atol=0.01) 
@test isapprox(LBO.n_Z_principal(0.5321u"µm", 293u"K"), 1.62122, atol=0.01) 
@test isapprox(LBO.n_Z_principal(1.1000u"µm", 293u"K"), 1.60449, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(LBO, 0.5321u"µm"), 109.2u"°"/2, atol=ustrip(u"rad", 1u"°")) 

# Test sampled phase-matches
pm1 = find_nearest_pm_along_theta_phi(90u"°", 11.3u"°", (:hi, :hi, :lo), LBO; lambda_r1=1.0642u"µm", lambda_b=0.5321u"µm", temp=293u"K")
@test all(pm1.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm1.pm_type[1].principal_plane == :XY
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.phi_pm, 11.3u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(abs(pm1.eff_data.d_eff), 0.831u"pm/V", rtol=0.1) # From SNLO

pm2 = find_nearest_pm_along_theta_phi(86.74u"°", 0u"°", (:hi, :hi, :lo), LBO; lambda_r1=1.3414u"µm", lambda_b=1.3414u"µm" / 2, temp=293u"K")
@test all(pm2.pm_type[1].o_or_e_rrb .== (:e, :e, :o))
@test pm2.pm_type[1].principal_plane == :XZ
@test isnothing(pm2.pm_type[2])
@test isapprox(pm2.theta_pm, 86.74u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(abs(pm2.eff_data.d_eff), 0.819u"pm/V", rtol=0.1) # From SNLO

pm3 = find_nearest_pm_along_theta_phi(90u"°", 21.11u"°", (:hi, :hi, :lo), LBO; lambda_r1=1.3188u"µm", lambda_b=0.4396u"µm", temp=293u"K")
@test all(pm3.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm3.pm_type[1].principal_plane == :XY
@test isnothing(pm3.pm_type[2])
@test isapprox(pm3.phi_pm, 21.11u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm3.walkoff_angle_rrb[1], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[3], 0.705u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(abs(pm3.eff_data.d_eff), 0.811u"pm/V", rtol=0.1) # From SNLO

pm4 = find_nearest_pm_along_theta_phi(5.1u"°", 0.0u"°", (:lo, :hi, :lo), LBO; lambda_r1=1.3188u"µm", lambda_b=0.6594u"µm", temp=293u"K")
@test all(pm4.pm_type[1].o_or_e_rrb .== (:e, :o, :e))
@test pm4.pm_type[1].principal_plane == :XZ
@test isnothing(pm4.pm_type[2])
@test isapprox(pm4.theta_pm, 5.1u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm4.walkoff_angle_rrb[1], 0.248u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_rrb[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_rrb[3], 0.262u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(abs(pm4.eff_data.d_eff), 0.642u"pm/V", rtol=0.1) # From SNLO

pm5 = find_nearest_pm_along_theta_phi(42.19u"°", 90.0u"°", (:lo, :hi, :lo), LBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test all(pm5.pm_type[1].o_or_e_rrb .== (:o, :e, :o))
@test pm5.pm_type[1].principal_plane == :YZ
@test isnothing(pm5.pm_type[2])
@test isapprox(pm5.theta_pm, 42.19u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm5.walkoff_angle_rrb[1], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm5.walkoff_angle_rrb[2], 0.533u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm5.walkoff_angle_rrb[3], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(abs(pm5.eff_data.d_eff), 0.536u"pm/V", rtol=0.1) # From SNLO

pm6 = find_all_pms_along_dimension((:hi, :hi, :lo), LBO; lambda_r1_fixed=1064u"nm", lambda_r2_fixed=1064u"nm", temp_min=273.15u"K", temp_max=600u"K", principal_axis=:X)[1]
@test isapprox(pm6.temp, 149u"°C", atol=5u"K")
@test isapprox(pm6.theta_pm, 90.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(pm6.phi_pm, 0.0u"°", atol=ustrip(u"rad", 0.1u"°"))
@test isapprox(abs(pm6.eff_data.d_eff), 0.85u"pm/V", rtol=0.1) # From SNLO

# Test sampled noncritical phase-matches
# TODO