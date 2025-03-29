using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_bbo()) <: UnidirectionalCrystal

@test BBO.n_o_principal(default_lambda(BBO), default_temp(BBO)) > 1
@test BBO.n_e_principal(default_lambda(BBO), default_temp(BBO)) > 1

# Sample refractive indices
@test isapprox(BBO.n_o_principal(0.40466u"µm", 293u"K"), 1.69267, atol=0.01) 
@test isapprox(BBO.n_o_principal(0.57907u"µm", 293u"K"), 1.67131, atol=0.01) 
@test isapprox(BBO.n_o_principal(1.01400u"µm", 293u"K"), 1.65608, atol=0.01) 

@test isapprox(BBO.n_e_principal(0.40466u"µm", 293u"K"), 1.56796, atol=0.01) 
@test isapprox(BBO.n_e_principal(0.57907u"µm", 293u"K"), 1.55298, atol=0.01) 
@test isapprox(BBO.n_e_principal(1.01400u"µm", 293u"K"), 1.54333, atol=0.01) 


# Test optical axis
@test isapprox(optical_axis_angle(BBO, 0.5321u"µm"), 0.0u"°", atol=ustrip(u"rad", 1u"°")) 

# Test sampled phasematches (from: Handbook of Nonlinear Crystals)
pm1 = find_nearest_pm_along_theta_phi(48u"°", 0.0u"°", [:o, :o, :e], BBO; lambda_r1=0.5321u"µm", lambda_b=0.26605u"µm", temp=293u"K")
@test isapprox(pm1.theta_pm, 47.62u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm1.walkoff_angle_r1_r2_b[1], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_r1_r2_b[2], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_r1_r2_b[3], 85.3u"mrad", atol=ustrip(u"rad", 5u"°")) 
# @test isapprox(abs(pm1.d_eff), 1.75u"pm/V", rtol=0.1) # From SNLO

pm2 = find_nearest_pm_along_theta_phi(38.39u"°", 0.0u"°", [:e, :o, :e], BBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test isapprox(pm2.theta_pm, 38.39u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm2.theta_L_bw, 0.02u"°" * 1u"cm", atol=0.1u"mrad * cm") 
@test isapprox(pm2.temp_L_bw, 16.37u"K" * 1u"cm", atol=2u"K * cm") # From SNLO
# @test isapprox(abs(pm2.d_eff), 1.29u"pm/V", rtol=0.1) # From SNLO

pm3 = find_nearest_pm_along_theta_phi(45.50u"°", 0.0u"°", [:o, :e, :e], BBO; lambda_r1=1.3188u"µm", lambda_b=0.4396u"µm", temp=293u"K")
@test isapprox(pm3.theta_pm, 45.50u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm3.walkoff_angle_r1_r2_b[1], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm3.walkoff_angle_r1_r2_b[2], 4.164u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm3.walkoff_angle_r1_r2_b[3], 4.312u"°", atol=ustrip(u"rad", 1u"°")) 
# @test isapprox(abs(pm3.d_eff), 0.909u"pm/V", rtol=0.1) # From SNLO

pm4 = find_nearest_pm_along_theta_phi(58.4u"°", 0.0u"°", [:o, :e, :e], BBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test isapprox(pm4.theta_pm, 58.4u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm4.theta_L_bw, 0.05u"°" * 1u"cm", atol=0.1u"mrad * cm")
@test isapprox(pm4.temp_L_bw, 18.43u"K" * 1u"cm", atol=2u"K * cm") # From SNLO
# @test isapprox(abs(pm4.d_eff), 0.478u"pm/V", rtol=0.1) # From SNLO
