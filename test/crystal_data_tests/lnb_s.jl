using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_lnb_s()) <: UnidirectionalCrystal
@test isa(crystal_system(LNB_S), String)

@test LNB_S.n_o_principal(default_lambda(LNB_S), default_temp(LNB_S)) > 1
@test LNB_S.n_e_principal(default_lambda(LNB_S), default_temp(LNB_S)) > 1

# Sample refractive indices (from: Handbook of Nonlinear Crystals)
@test isapprox(LNB_S.n_o_principal(0.42u"µm", 293u"K"), 2.4089, atol=0.01) 
@test isapprox(LNB_S.n_o_principal(1.60u"µm", 293u"K"), 2.2113, atol=0.01) 
@test isapprox(LNB_S.n_o_principal(3.0u"µm", 293u"K"), 2.1625, atol=0.01) 
@test isapprox(LNB_S.n_o_principal(4.0u"µm", 293u"K"), 2.1155, atol=0.01) 

@test isapprox(LNB_S.n_e_principal(0.42u"µm", 293u"K"), 2.3025, atol=0.01) 
@test isapprox(LNB_S.n_e_principal(1.60u"µm", 293u"K"), 2.1361, atol=0.01) 
@test isapprox(LNB_S.n_e_principal(3.0u"µm", 293u"K"), 2.0945, atol=0.01) 
@test isapprox(LNB_S.n_e_principal(4.0u"µm", 293u"K"), 2.0553, atol=0.01) 

# Test optical axis
@test isapprox(optical_axis_angle(LNB_S, 0.5321u"µm"), 0.0u"°", atol=ustrip(u"rad", 1u"°")) 

# Test sampled phasematches
pm1 = find_nearest_pm_along_theta_phi(71.8u"°", 30.0u"°", (:hi, :hi, :lo), LNB_S; lambda_r1=1.118u"µm", lambda_b=0.559u"µm", temp=293u"K")
@test all(pm1.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm1.pm_type[1].principal_plane == :UD
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.theta_pm, 71.8u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(pm1.phi_pm, 30.0u"°", atol=ustrip(u"rad", 1u"°"))
@test isapprox(abs(pm1.eff_data.d_eff), 4.94u"pm/V", rtol=0.1) # From SNLO

# Test sampled noncritical phasematches 
# TODO