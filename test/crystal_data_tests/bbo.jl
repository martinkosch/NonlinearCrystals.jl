using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_bbo()) <: UnidirectionalCrystal
@test isa(crystal_system(BBO), String)

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

# Test sampled phase-matches 
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#nikogosyan91beta
pm1 = find_nearest_pm_along_theta_phi(48u"°", 90.0u"°", (:o, :o, :e), BBO; lambda_r1=0.5321u"µm", lambda_b=0.26605u"µm", temp=293u"K")
@test isapprox(pm1.theta_pm, 47.62u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[1], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[2], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[3], 85.3u"mrad", atol=ustrip(u"rad", 5u"°")) 
@test isapprox(abs(pm1.eff_data.d_eff), 1.62u"pm/V", rtol=0.2) # From Nikogosyan

pm2 = find_nearest_pm_along_theta_phi(38.39u"°", 60.0u"°", (:e, :o, :e), BBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test isapprox(pm2.theta_pm, 38.39u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm2.bw_data.theta_L_bw, 0.02u"°" * 1u"cm", atol=0.1u"mrad * cm") 
@test isapprox(pm2.bw_data.temp_L_bw, 16.37u"K" * 1u"cm", atol=2u"K * cm") # From SNLO
@test isapprox(abs(pm2.eff_data.d_eff), 1.37u"pm/V", rtol=0.2) # From Nikogosyan

pm3 = find_nearest_pm_along_theta_phi(45.50u"°", 0.0u"°", (:o, :e, :e), BBO; lambda_r1=1.3188u"µm", lambda_b=0.4396u"µm", temp=293u"K")
@test isapprox(pm3.theta_pm, 45.50u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[1], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[2], 4.164u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[3], 4.312u"°", atol=ustrip(u"rad", 1u"°")) 

pm4 = find_nearest_pm_along_theta_phi(58.4u"°", 0.0u"°", (:o, :e, :e), BBO; lambda_r1=1.0642u"µm", lambda_b=0.35473u"µm", temp=293u"K")
@test isapprox(pm4.theta_pm, 58.4u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm4.bw_data.theta_L_bw, 0.05u"°" * 1u"cm", atol=0.1u"mrad * cm")
@test isapprox(pm4.bw_data.temp_L_bw, 18.43u"K" * 1u"cm", atol=2u"K * cm") # From SNLO
@test isapprox(abs(pm4.eff_data.d_eff), 0.56u"pm/V", rtol=0.3) # From SNLO

# TODO: This is not working yet, SNLO seems to use non-textbook Miller scaling?
# pm5 = find_nearest_pm_along_theta_phi(69.5u"°", 30.0u"°", (:o, :o, :e), BBO; lambda_r1=0.5321u"µm", lambda_b=0.2128u"µm", temp=293u"K")
# @test isapprox(pm5.theta_pm, 69.5u"°", atol=ustrip(u"rad", 2u"°")) 
# @test isapprox(abs(pm5.eff_data.d_eff), 0.93u"pm/V", rtol=0.2) # From Nikogosyan; SNLO: 1.09u"pm/V"
