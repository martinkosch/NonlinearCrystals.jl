using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_lnb_m()) <: UnidirectionalCrystal
@test isa(crystal_system(LNB_M), String)

@test LNB_M.n_o_principal(default_lambda(LNB_M), default_temp(LNB_M)) > 1
@test LNB_M.n_e_principal(default_lambda(LNB_M), default_temp(LNB_M)) > 1

# Sample refractive indices
# TODO

# Test optical axis
@test isapprox(optical_axis_angle(LNB_M, 0.5321u"µm"), 0.0u"°", atol=ustrip(u"rad", 1u"°")) 

# Test sampled phase-matches
# TODO

# Test sampled noncritical phase-matches
# TODO