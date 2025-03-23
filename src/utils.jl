export all_nonlinearcrystals, all_unidirectionalcrystals, all_bidirectionalcrystals

function subtypes_in_module(mod::Module, T; recurse=false)
    found = T[]
    for name in names(mod)
        binding = getfield(mod, name)
        if binding isa T
            push!(found, binding)
        elseif recurse && binding isa Module
            append!(found, subtypes_in_module(binding, T; recurse=true))
        end
    end
    return found
end

function all_nonlinearcrystals()
    return subtypes_in_module(NonlinearCrystals, NonlinearCrystal)
end

function all_unidirectionalcrystals()
    return subtypes_in_module(NonlinearCrystals, UnidirectionalCrystal)
end

function all_bidirectionalcrystals()
    return subtypes_in_module(NonlinearCrystals, BidirectionalCrystal)
end

function draw_ellipsoid!(
    ax::Axis3,
    center::AbstractVector=[0.0, 0.0, 0.0],
    principal_axes::AbstractVector=[1.0, 1.0, 1.0];
    kwargs...,
)
    s = Sphere(Point3f(center...), 1f0)

    meshplot = mesh!(ax, s, alpha=0.3, transparency=true)
    scale!(meshplot, principal_axes...)

    return ax
end

function vector_to_angles(v::AbstractVector{<:Number})
    x, y, z = v
    r = sqrt(x^2 + y^2 + z^2)
    θ = Float64(acos(z / r)) |> u"rad" # Polar angle (0 to π)
    ϕ = atan(y, x)  # Azimuthal angle (−π to π)
    ϕ = mod(ϕ, 2π) |> u"rad" # Normalize φ to [0, 2π)
    return θ, ϕ
end

function angles_to_vector(theta::Angle, phi::Angle)
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return @SVector [x, y, z]
end

function axes_to_θ_ϕ(axes::Union{AbstractVector{Symbol}, Symbol})
    axes = (axes isa Symbol) ? [axes] : axes
    θ_ϕ_dirs = map(axes) do a
        if a == :X
            return (90u"°", 0u"°")
        elseif a==:Y
            return (90u"°", 90u"°")
        elseif a==:Z
            return (0u"°", 0u"°")
        else
            @error "Axis '$(a)' unknown. Valid options are ':X', ':Y', and ':Z'."
        end
    end
    return θ_ϕ_dirs
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

function shift_lambda_with_freq(lambda::Length, delta_omega::Frequency)
    omega = c_0 ./ (2π * lambda)
    lambda_shifted = c_0 / (2π * (omega + delta_omega))
    return lambda_shifted
end