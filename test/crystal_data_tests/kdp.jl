using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_kdp()) <: UnidirectionalCrystal
@test isa(crystal_system(KDP), String)

@test KDP.n_o_principal(default_lambda(KDP), default_temp(KDP)) > 1
@test KDP.n_e_principal(default_lambda(KDP), default_temp(KDP)) > 1

# Sample refractive indices
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
@test isapprox(KDP.n_o_principal(0.2138560u"µm", 298u"K"), 1.60177, atol=0.01) 
@test isapprox(KDP.n_o_principal(0.5460740u"µm", 298u"K"), 1.51152, atol=0.01) 
@test isapprox(KDP.n_o_principal(1.1522760u"µm", 298u"K"), 1.49135, atol=0.01) 

@test isapprox(KDP.n_e_principal(0.2138560u"µm", 298u"K"), 1.54615, atol=0.01) 
@test isapprox(KDP.n_e_principal(0.5460740u"µm", 298u"K"), 1.46982, atol=0.01) 
@test isapprox(KDP.n_e_principal(1.1522760u"µm", 298u"K"), 1.45893, atol=0.01) 


# Test temperature derivatives
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
@test isapprox(dn_dtemp(KDP.n_o_principal, 0.405u"µm", 293u"K"), -3.27e-5u"K^-1", atol=1.0e-5u"K^-1")
@test isapprox(dn_dtemp(KDP.n_o_principal, 0.633u"µm", 293u"K"), -3.94e-5u"K^-1", atol=1.0e-5u"K^-1")

@test isapprox(dn_dtemp(KDP.n_e_principal, 0.405u"µm", 293u"K"), -3.15e-5u"K^-1", atol=1.0e-5u"K^-1")
@test isapprox(dn_dtemp(KDP.n_e_principal, 0.633u"µm", 293u"K"), -2.54e-5u"K^-1", atol=1.0e-5u"K^-1")

# Test sampled phase-matches 
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
pm1 = find_nearest_pm_along_theta_phi(40u"°", 90.0u"°", (:o, :o, :e), KDP; lambda_r1=1.06u"µm", lambda_r2=1.06u"µm", temp=293u"K")
@test isapprox(pm1.theta_pm, 41u"°", atol=ustrip(u"rad", 2u"°")) 

pm2 = find_nearest_pm_along_theta_phi(55u"°", 90.0u"°", (:e, :o, :e), KDP; lambda_r1=1.0642u"µm", lambda_r2=0.5321u"µm", temp=293u"K")
@test isapprox(pm2.theta_pm, 58.3u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm2.walkoff_angle_rrb[1], 1.149u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm2.walkoff_angle_rrb[2], 0.0u"°", atol=ustrip(u"rad", 1u"°")) 
@test isapprox(pm2.walkoff_angle_rrb[3], 1.404u"°", atol=ustrip(u"rad", 5u"°")) 
