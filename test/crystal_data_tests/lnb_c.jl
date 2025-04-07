using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_lnb_c()) <: UnidirectionalCrystal
@test isa(crystal_system(LNB_C), String)

@test LNB_C.n_o_principal(default_lambda(LNB_C), default_temp(LNB_C)) > 1
@test LNB_C.n_e_principal(default_lambda(LNB_C), default_temp(LNB_C)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(LNB_C.n_o_principal(0.43584u"µm", 293u"K"), 2.39276, atol=0.01) 
@test isapprox(LNB_C.n_o_principal(0.63282u"µm", 293u"K"), 2.28647, atol=0.01) 
@test isapprox(LNB_C.n_o_principal(3.3913u"µm", 293u"K"), 2.1415, atol=0.01) 

@test isapprox(LNB_C.n_e_principal(0.43584u"µm", 293u"K"), 2.29278, atol=0.01) 
@test isapprox(LNB_C.n_e_principal(0.63282u"µm", 293u"K"), 2.20240, atol=0.01) 
@test isapprox(LNB_C.n_e_principal(3.3913u"µm", 293u"K"), 2.0822, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(LNB_C, 0.5321u"µm"), 0.0u"°", atol=ustrip(u"rad", 1u"°")) 

# Test sampled phasematches
pm1 = find_nearest_pm_along_theta_phi(45u"°", 30.0u"°", (:hi, :hi, :lo), LNB_C; lambda_r1=2.12u"µm", lambda_b=1.06u"µm", temp=293u"K")
@test all(pm1.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm1.pm_type[1].principal_plane == :UD
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.theta_pm, 45.25u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm1.phi_pm, 30.0u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm1.walkoff_angle_rrb[1], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[2], 0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[3], 1.988u"°", atol=ustrip(u"rad", 0.2u"°")) 
@test isapprox(abs(pm1.eff_data.d_eff), 4.03u"pm/V", rtol=0.2) # From SNLO

# Test sampled noncritical phasematches 
# TODO