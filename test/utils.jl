using NonlinearCrystals
using Test
using Unitful
using LinearAlgebra

vectors = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0]
]

for v in vectors
    θ, ϕ = NonlinearCrystals.vector_to_angles(v)
    v_reconstructed = NonlinearCrystals.angles_to_vector(θ, ϕ)
    v_norm = normalize(v)
    @test isapprox(v_reconstructed, v_norm, atol=1e-10)
end


θϕ_X = NonlinearCrystals.axes_to_θ_ϕ(:X)[1]
@test isapprox(θϕ_X[1], 90u"°", atol=1e-8u"°")
@test isapprox(θϕ_X[2], 0u"°", atol=1e-8u"°")

θϕ_Y = NonlinearCrystals.axes_to_θ_ϕ(:Y)[1]
@test isapprox(θϕ_Y[1], 90u"°", atol=1e-8u"°")
@test isapprox(θϕ_Y[2], 90u"°", atol=1e-8u"°")

θϕ_Z = NonlinearCrystals.axes_to_θ_ϕ(:Z)[1]
@test isapprox(θϕ_Z[1], 0u"°", atol=1e-8u"°")
@test isapprox(θϕ_Z[2], 0u"°", atol=1e-8u"°")


perms = NonlinearCrystals.bool_permutations(:A, :B, 2)
expected = [
    (:A, :A),
    (:B, :A),
    (:A, :B),
    (:B, :B)
]
@test length(perms) == 4
@test all(p -> p in expected, perms)

## Custom Eigen tests
# --- Diagonal matrix (distinct)
A = [2.0 0.0; 0.0 3.0]
vals, vecs = NonlinearCrystals.eigen_2d(A)
@test sort(vals) ≈ [2.0, 3.0]
@test isapprox(dot(vecs[:,1], vecs[:,2]), 0.0; atol=1e-8)

# --- Identity matrix (degenerate eigenvalues)
A = [1.0 0.0; 0.0 1.0]
vals, vecs = NonlinearCrystals.eigen_2d(A)
@test sort(vals) ≈ [1.0, 1.0]
@test isapprox(norm(vecs[:,1]), 1.0; atol=1e-8)
@test isapprox(norm(vecs[:,2]), 1.0; atol=1e-8)

# --- Symmetric off-diagonal
A = [2.0 1.0; 1.0 2.0]
vals, vecs = NonlinearCrystals.eigen_2d(A)
λ_expected = [1.0, 3.0]
@test sort(vals) ≈ λ_expected
@test isapprox(dot(vecs[:,1], vecs[:,2]), 0.0; atol=1e-8)

# --- Complex eigenvalues but the function is only suited for symmetric matrices 
A = [0.0 -1.0; 1.0 0.0]
@test_throws AssertionError NonlinearCrystals.eigenvals_2d(A)

# --- Negative eigenvalues
A = [-2.0 0.0; 0.0 -3.0]
vals, vecs = NonlinearCrystals.eigen_2d(A)
@test sort(vals) ≈ [-3.0, -2.0]
@test isapprox(dot(vecs[:,1], vecs[:,2]), 0.0; atol=1e-8)


λ0 = 1.0u"μm"
Δω = 1e12u"Hz"
λ_shifted = NonlinearCrystals.shift_lambda_with_freq(λ0, Δω)
@test λ_shifted < λ0  # Positive Δω reduces wavelength