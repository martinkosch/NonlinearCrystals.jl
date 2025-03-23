using NonlinearCrystals
using Test
using Unitful
using GLMakie

# Define a simple SellmeierFunction for testing
n_fun = (λ, T) -> 1.5 + 0.01 * (λ / 1u"µm")  # Arbitrary linear function
dn_dtemp_fun = (λ, T) -> 1e-4u"K^-1"            # Constant derivative
sri = SellmeierFunction(n_fun, (0.5u"µm", 2.0u"µm"); dn_dtemp_fun=dn_dtemp_fun, temp_range=(200u"K", 400u"K"))

# Test refractive_index basic
λ_test = 1.0u"µm"
T_test = 300u"K"
n_val = refractive_index(sri, λ_test, T_test)
@test isapprox(n_val, 1.51 + 1e-4 * (300 - 293.15); atol=1e-6)

# Test lambda range validity
@test is_lambda_valid(1.0u"µm", sri) == true
@test is_lambda_valid(2.5u"µm", sri) == false

# Test dn/dλ
dn = NonlinearCrystals.dn_dλ(sri, 1.0u"µm", T_test)
@test dn ≈ 0.01u"µm^-1" atol=1e-4u"m^-1"

# Test d²n/dλ² (should be ~0 because n_fun linear)
d2n = NonlinearCrystals.d2n_dλ2(sri, 1.0u"µm", T_test)
@test abs(d2n) < 1e-8u"m^-2"

# Test β0
β0_val = β0(sri, λ_test, T_test)
@test β0_val isa Unitful.Quantity

# Test β1 (group delay)
β1_val = β1(sri, λ_test, T_test)
@test β1_val isa Unitful.Quantity

# Test β2
β2_val = β2(sri, λ_test, T_test)
@test β2_val isa Unitful.Quantity

# Test β3 existence
β3_val = β3(sri, λ_test, T_test)
@test β3_val isa Unitful.Quantity

# Test group index and velocities
g_idx = group_index(sri, λ_test, T_test)
@test g_idx isa Real

v_p = phase_velocity(sri, λ_test, T_test)
v_g = group_velocity(sri, λ_test, T_test)
@test v_p > 0u"m/s"
@test v_g > 0u"m/s"

# Plot function (only test that it runs without error)
f = plot_refractiveindex(sri; temp=[293.15u"K", 310u"K"])
@test f isa Figure
