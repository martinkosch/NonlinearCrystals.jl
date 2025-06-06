using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_gs()) <: UnidirectionalCrystal
@test isa(crystal_system(GS), String)

@test GS.n_o_principal(default_lambda(GS), default_temp(GS)) > 1
@test GS.n_e_principal(default_lambda(GS), default_temp(GS)) > 1

# Sample refractive indices
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#vodopyanov1995
# Data read from diagrams; no tabulated data available in this reference
@test isapprox(GS.n_o_principal(0.6328u"µm", 293u"K"), 2.95, atol=0.02) 
@test isapprox(GS.n_o_principal(1.1523u"µm", 293u"K"), 2.78, atol=0.02) 
@test isapprox(GS.n_o_principal(3.3913u"µm", 293u"K"), 2.74, atol=0.02) 

@test isapprox(GS.n_e_principal(0.6328u"µm", 293u"K"), 2.7, atol=0.02) 
@test isapprox(GS.n_e_principal(1.1523u"µm", 293u"K"), 2.44, atol=0.02) 
@test isapprox(GS.n_e_principal(3.3913u"µm", 293u"K"), 2.4, atol=0.02) 

# Test sampled phase-matches 
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
pm1 = find_nearest_pm_along_theta_phi(18u"°", 30u"°", (:o, :o, :e), GS; lambda_r1=2.36u"µm", lambda_r2=2.36u"µm", temp=293u"K") 
@show pm1
@test all(pm1.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm1.pm_type[1].principal_plane == :UD
@test isnothing(pm1.pm_type[2])
@test isapprox(pm1.theta_pm, 18.7u"°", atol=ustrip(u"rad", 2u"°"))
@test isapprox(pm1.walkoff_angle_rrb[1], 0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm1.walkoff_angle_rrb[2], 0u"°", atol=ustrip(u"rad", 1u"°")) 

pm2 = find_nearest_pm_along_theta_phi(13u"°", 30u"°", (:o, :o, :e), GS; lambda_r1=10.6u"µm", lambda_r2=10.6u"µm", temp=293u"K") 
@test all(pm2.pm_type[1].o_or_e_rrb .== (:o, :o, :e))
@test pm2.pm_type[1].principal_plane == :UD
@test isnothing(pm2.pm_type[2])
@test isapprox(pm2.theta_pm, 14.89u"°", atol=ustrip(u"rad", 2u"°"))
@test isapprox(pm2.walkoff_angle_rrb[1], 0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm2.walkoff_angle_rrb[2], 0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm2.walkoff_angle_rrb[3], 4.059u"°", atol=ustrip(u"rad", 1u"°")) 

pm3 = find_nearest_pm_along_theta_phi(13u"°", 30u"°", (:o, :e, :e), GS; lambda_r1=9.6u"µm", lambda_r2=4.8u"µm", temp=293u"K") 
@test all(pm3.pm_type[1].o_or_e_rrb .== (:o, :e, :e))
@test pm3.pm_type[1].principal_plane == :UD
@test isnothing(pm3.pm_type[2])
@test isapprox(pm3.theta_pm, 18.57u"°", atol=ustrip(u"rad", 2u"°"))
@test isapprox(pm3.walkoff_angle_rrb[1], 0u"°", atol=ustrip(u"rad", 0.1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[2], 4.881u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm3.walkoff_angle_rrb[3], 4.876u"°", atol=ustrip(u"rad", 1u"°")) 
