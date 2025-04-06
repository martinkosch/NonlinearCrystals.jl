using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_bibo()) <: BidirectionalCrystal
@test isa(crystal_system(BIBO), String)

@test BIBO.n_X_principal(default_lambda(BIBO), default_temp(BIBO)) > 1
@test BIBO.n_Y_principal(default_lambda(BIBO), default_temp(BIBO)) > 1
@test BIBO.n_Z_principal(default_lambda(BIBO), default_temp(BIBO)) > 1

# Sample refractive indices (K. Miyata, N. Umemura, and K. Kato, Phase-matched pure χ(3) third-harmonic generation in noncentrosymmetric BiB3O6 , Opt. Lett. 34, 500–502 (2009).)
@test isapprox(BIBO.n_X_principal(0.6328u"µm", 293u"K"), 1.77668, atol=0.01)
@test isapprox(BIBO.n_X_principal(1.064u"µm", 293u"K"), 1.75752, atol=0.01)

@test isapprox(BIBO.n_Y_principal(0.6328u"µm", 293u"K"), 1.80641, atol=0.01)
@test isapprox(BIBO.n_Y_principal(1.064u"µm", 293u"K"), 1.78400, atol=0.01)

@test isapprox(BIBO.n_Z_principal(0.6328u"µm", 293u"K"), 1.94582, atol=0.01)
@test isapprox(BIBO.n_Z_principal(1.064u"µm", 293u"K"), 1.91711, atol=0.01)

# Test temperature derivatives
@test isapprox(dn_dtemp(BIBO.n_X_principal, 0.6328u"µm", 293u"K"), 8.6e-6u"K^-1", atol=1.0e-6u"K^-1")
@test isapprox(dn_dtemp(BIBO.n_X_principal, 1.064u"µm", 293u"K"), 3.5e-6u"K^-1", atol=1.0e-6u"K^-1")

@test isapprox(dn_dtemp(BIBO.n_Y_principal, 0.6328u"µm", 293u"K"), -3.7e-6u"K^-1", atol=1.0e-6u"K^-1")
@test isapprox(dn_dtemp(BIBO.n_Y_principal, 1.064u"µm", 293u"K"), -5.6e-6u"K^-1", atol=1.0e-6u"K^-1")

@test isapprox(dn_dtemp(BIBO.n_Z_principal, 0.6328u"µm", 293u"K"), -6.3e-6u"K^-1", atol=1.0e-6u"K^-1")
@test isapprox(dn_dtemp(BIBO.n_Z_principal, 1.064u"µm", 293u"K"), -6.8e-6u"K^-1", atol=1.0e-6u"K^-1")


# Test sampled phasematches
pm1 = find_nearest_pm_along_theta_phi(150.0u"°", 90.0u"°", (:hi, :hi, :lo), BIBO; lambda_r1=0.8u"µm", lambda_b=0.40u"µm", temp=293u"K")
@test isapprox(pm1.theta_pm, 151.2u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm1.phi_pm, 90.0u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(abs(pm1.eff_data.d_eff), 3.59u"pm/V", rtol=0.2) # From SNLO 

pm2 = find_nearest_pm_along_theta_phi(51.5u"°", 0.0u"°", (:hi, :hi, :lo), BIBO; lambda_r1=0.8u"µm", lambda_b=0.40u"µm", temp=293u"K")
@test isapprox(pm2.theta_pm, 51.5u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm2.phi_pm, 0.0u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(abs(pm2.eff_data.d_eff), 1.81u"pm/V", rtol=0.2) # From SNLO 


pm3 = find_nearest_pm_along_theta_phi(11.1u"°", 0.0u"°", (:hi, :hi, :lo), BIBO; lambda_r1=1.6u"µm", lambda_b=0.80u"µm", temp=293u"K")
@test isapprox(pm3.theta_pm, 11.1u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(pm3.phi_pm, 0.0u"°", atol=ustrip(u"rad", 2u"°")) 
@test isapprox(abs(pm3.eff_data.d_eff), 2.69u"pm/V", rtol=0.2) # From SNLO 

# Test sampled noncritical phasematches
# TODO