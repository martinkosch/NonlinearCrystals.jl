using NonlinearCrystals
using Test
using Unitful
using GLMakie
import PhysicalConstants.CODATA2022: c_0

# Define a simple SellmeierFunction for testing
temp_ref = 293.15u"K"
n_fun = (λ, T) -> 1.5 + 0.01 * (λ / 1u"µm") + 1e-4u"K^-1" * (T - temp_ref) # Arbitrary linear function
sri = SellmeierFunction(n_fun, (0.5u"µm", 2.0u"µm"); temp_ref, temp_range=(200u"K", 400u"K"))

# Test refractive_index basic
λ_test = 1.0u"µm"
T_test = 300u"K"
n_val = refractive_index(sri, λ_test, T_test)
@test isapprox(n_val, 1.51 + 1e-4u"K^-1" * (T_test - temp_ref); atol=1e-6)

# Test lambda range validity
@test is_lambda_valid(1.0u"µm", sri) == true
@test is_lambda_valid(2.5u"µm", sri) == false

# Test dn/dλ
dn = NonlinearCrystals.dn_dλ(sri, 1.0u"µm", T_test)
@test dn ≈ 0.01u"µm^-1" atol = 1e-4u"m^-1"

# Test d²n/dλ² (should be ~0 because n_fun linear)
d2n = NonlinearCrystals.d2n_dλ2(sri, 1.0u"µm", T_test)
@test abs(d2n) < 1e-8u"m^-2"

# Test β0
β0_val = β0(sri, λ_test, T_test)
@test isapprox(uconvert(u"m^-1", β0_val), 2π / (1.0u"µm") * (1.5 + 0.01), rtol=0.01)  # Expected k ≈ 2π/λ × n

# Test β1 (group delay)
β1_val = β1(sri, λ_test, T_test)
@test isapprox(uconvert(u"s/m", β1_val), 1.5 / c_0, rtol=0.01)

# Test β2
β2_val = β2(sri, λ_test, T_test)
@test isapprox(uconvert(u"s^2/m", β2_val), 0.0u"s^2/m", atol=0.01u"fs^2/mm")

# Test β3 existence
β3_val = β3(sri, λ_test, T_test)
@test isapprox(uconvert(u"s^3/m", β3_val), 0.0u"s^3/m", atol=0.01u"fs^3/mm")

# Test group velocity index and velocities
g_idx = group_index(sri, λ_test, T_test)
@test unit(g_idx) == Unitful.NoUnits
@test g_idx > 0.0  # physically meaningful
@test g_idx < 100  # physically meaningful

v_p = phase_velocity(sri, λ_test, T_test)
v_g = group_velocity(sri, λ_test, T_test)
@test v_p > 0u"m/s"
@test v_p < c_0
@test v_g > 0u"m/s"

# Plot function (only test that it runs without error)
f = plot_refractiveindex(sri; temp=[293.15u"K", 310u"K"])
@test f isa Figure
