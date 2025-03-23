using NonlinearCrystals
using Test
using Unitful

@test LBO.n_x_principal(default_lambda(LBO), default_temp(LBO)) > 1
@test LBO.n_y_principal(default_lambda(LBO), default_temp(LBO)) > 1
@test LBO.n_z_principal(default_lambda(LBO), default_temp(LBO)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(LBO.n_x_principal(0.2537u"µm", 293u"K"), 1.6335, atol=0.01) 
@test isapprox(LBO.n_x_principal(0.5321u"µm", 293u"K"), 1.57868, atol=0.01) 
@test isapprox(LBO.n_x_principal(1.1000u"µm", 293u"K"), 1.56432, atol=0.01) 

@test isapprox(LBO.n_y_principal(0.2537u"µm", 293u"K"), 1.6582, atol=0.01) 
@test isapprox(LBO.n_y_principal(0.5321u"µm", 293u"K"), 1.60642, atol=0.01) 
@test isapprox(LBO.n_y_principal(1.1000u"µm", 293u"K"), 1.59005, atol=0.01) 

@test isapprox(LBO.n_z_principal(0.2537u"µm", 293u"K"), 1.6792, atol=0.01) 
@test isapprox(LBO.n_z_principal(0.5321u"µm", 293u"K"), 1.62122, atol=0.01) 
@test isapprox(LBO.n_z_principal(1.1000u"µm", 293u"K"), 1.60449, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(LBO, 0.5321u"µm"), 109.2u"°"/2, atol=ustrip(u"rad", 1u"°")) 

# Test sampled phasematches
pm1 = find_nearest_pm(90u"°", 11.3u"°", [:hi, :hi, :lo], LBO; lambda_r1=1.0642u"µm", lambda_b=0.5321u"µm", temp=293u"K")
@test all(isnothing.(pm1.o_or_e_r1_r2_b)) || all(pm1.o_or_e_r1_r2_b .== [:o, :o, :e])
@test isapprox(pm1.phi_pm, 11.3u"°", atol=ustrip(u"rad", 1u"°"))

pm2 = find_nearest_pm(86.74u"°", 0u"°", [:hi, :hi, :lo], LBO; lambda_r1=1.3414u"µm", lambda_b=1.3414u"µm" / 2, temp=293u"K")
@test all(isnothing.(pm2.o_or_e_r1_r2_b)) || all(pm2.o_or_e_r1_r2_b .== [:e, :e, :o])
@test isapprox(pm2.theta_pm, 86.74u"°", atol=ustrip(u"rad", 1u"°"))

pm3 = find_nearest_pm(90u"°", 21.11u"°", [:hi, :hi, :lo], LBO; lambda_r1=1.3188u"µm", lambda_b=0.4396u"µm", temp=293u"K")
@test all(isnothing.(pm3.o_or_e_r1_r2_b)) || all(pm3.o_or_e_r1_r2_b .== [:o, :o, :e])
@test isapprox(pm3.phi_pm, 21.11u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm3.walkoff_angle_r1_r2_b[1], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm3.walkoff_angle_r1_r2_b[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm3.walkoff_angle_r1_r2_b[3], 0.705u"°", atol=ustrip(u"rad", 0.1u"°")) 

pm4 = find_nearest_pm(5.1u"°", 0.0u"°", [:lo, :hi, :lo], LBO; lambda_r1=1.3188u"µm", lambda_b=0.6594u"µm", temp=293u"K")
@test all(isnothing.(pm4.o_or_e_r1_r2_b)) || all(pm4.o_or_e_r1_r2_b .== [:e, :o, :e])
@test isapprox(pm4.theta_pm, 5.1u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm4.walkoff_angle_r1_r2_b[1], 0.248u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_r1_r2_b[2], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_r1_r2_b[3], 0.262u"°", atol=ustrip(u"rad", 0.1u"°")) 

pm4 = find_nearest_pm(42.19u"°", 0.0u"°", [:lo, :hi, :lo], LBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test all(isnothing.(pm4.o_or_e_r1_r2_b)) || all(pm4.o_or_e_r1_r2_b .== [:o, :e, :o])
@test isapprox(pm4.theta_pm, 42.19u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm4.walkoff_angle_r1_r2_b[1], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_r1_r2_b[2], 0.533u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm4.walkoff_angle_r1_r2_b[3], 0.0u"°", atol=ustrip(u"rad", 0.1u"°")) 
