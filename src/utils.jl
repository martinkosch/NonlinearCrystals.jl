export draw_ellipsoid!

function draw_ellipsoid!(
    ax,
    center::AbstractVector=[0.0, 0.0, 0.0],
    principal_axes::AbstractVector=[1.0, 1.0, 1.0];
    # rotation::AbstractVector=[0.0, 0.0, 0.0];
    kwargs...,
)
    s = Sphere(Point3f(center...), 1f0)

    meshplot = mesh!(ax, s, alpha=0.3, transparency=true)
    scale!(meshplot, principal_axes...)

    return ax
end

function vector_to_angles(v::AbstractVector{<:Real})
    x, y, z = v
    r = sqrt(x^2 + y^2 + z^2)
    θ = Float64(acos(z / r)) |> u"rad"   # Polar angle (0 to π)
    ϕ = atan(y, x)  # Azimuthal angle (−π to π)
    ϕ = mod(ϕ, 2π) |> u"rad"    # Normalize φ to [0, 2π)
    return θ, ϕ
end

function angles_to_vector(θ, ϕ)
    x = sin(θ) * cos(ϕ)
    y = sin(θ) * sin(ϕ)
    z = cos(θ)
    return @SVector [x, y, z]
end

function bool_permutations(val1, val2, n::Integer)
    result = []

    # Iterate over all possible 2^n combinations
    for i in 0:(2^n - 1)
        perm = Vector{typeof(val1)}(undef, n)
        for j in 1:n
            # Choose val1 or val2 based on bit pattern
            if (i >> (j - 1)) & 1 == 1
                perm[j] = val2
            else
                perm[j] = val1
            end
        end
        push!(result, perm)
    end

    return result
end

function eigen_2d(A::AbstractMatrix)
    (a, b, c, d) = A'
    discr = sqrt(max(0.0, ((a + d)/2)^2 - (a * d - b * c))) # Bound from below to avoid numerical issues
    eigvals = (a + d)/2 .+ [-discr, discr]

    eigvecs_orth1 = [[b, eigvals[1] - a], [eigvals[1] - d, c]]
    idx_ev1 = findmax(norm.(eigvecs_orth1))[2]
    eigvecs_orth2 = [[b, eigvals[2] - a], [eigvals[2] - d, c]]
    idx_ev2 = findmax(norm.(eigvecs_orth2))[2]
    eigvecs_orth = [normalize(eigvecs_orth1[idx_ev1]) normalize(eigvecs_orth2[idx_ev2])]
    return eigvals, eigvecs_orth
end