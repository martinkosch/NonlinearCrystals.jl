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

function auto_fmt(x; digits=3, sci_thresh=1e4, inf_thresh=1e10)
    if isnan(x)
        return "NaN"
    elseif isinf(x) || abs(x) ≥ inf_thresh
        return x > 0 ? "Inf" : "-Inf"
    end

    absx = abs(x)
    if absx ≥ sci_thresh || (absx > 0 && absx < 0.01)
        return @sprintf("%.*e", digits, x)
    else
        return @sprintf("%.*f", digits, x)
    end
end

function vec_str(v::AbstractVector{<:Number}; digits::Integer=3) 
    return "[" * join(round.(v; digits), ", ") * "]"
end

function draw_ellipsoid!(
    ax::Axis3,
    center::AbstractVector=[0.0, 0.0, 0.0],
    principal_axes::AbstractVector=[1.0, 1.0, 1.0];
    kwargs...,
)
    s = Sphere(Point3f(center...), 1.0f0)

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

function axes_to_θ_ϕ(axes::Union{AbstractVector{Symbol},Symbol})
    axes = (axes isa Symbol) ? [axes] : axes
    θ_ϕ_dirs = map(axes) do a
        if a == :X
            return (90u"°", 0u"°")
        elseif a == :Y
            return (90u"°", 90u"°")
        elseif a == :Z
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
    for i in 0:(2^n-1)
        perm = Vector{typeof(val1)}(undef, n)
        for j in 1:n
            # Choose val1 or val2 based on bit pattern
            if (i >> (j - 1)) & 1 == 1
                perm[j] = val2
            else
                perm[j] = val1
            end
        end
        push!(result, Tuple(perm))
    end

    return result
end

function eigen_2d(A::AbstractMatrix)
    (a, b, c, d) = A'
    discr = sqrt(max(0.0, ((a + d) / 2)^2 - (a * d - b * c))) # Bound from below to avoid numerical issues
    eigvals = (a + d) / 2 .+ [-discr, discr]

    eigvecs_orth1 = [(@SVector [b, eigvals[1] - a]), (@SVector [eigvals[1] - d, c])]
    idx_ev1 = findmax(norm.(eigvecs_orth1))[2]
    eigvecs_orth2 = [(@SVector [b, eigvals[2] - a]), (@SVector [eigvals[2] - d, c])]
    idx_ev2 = findmax(norm.(eigvecs_orth2))[2]
    eigvecs_orth = [normalize(eigvecs_orth1[idx_ev1]) normalize(eigvecs_orth2[idx_ev2])]
    return eigvals, eigvecs_orth
end

function shift_lambda_with_freq(lambda::Length, delta_omega::Frequency)
    omega = c_0 ./ (2π * lambda)
    lambda_shifted = c_0 / (2π * (omega + delta_omega))
    return lambda_shifted
end

function split_on_nan(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    @assert length(x) == length(y) "Vectors must be the same length"

    segments_x = Vector{T}[]
    segments_y = Vector{T}[]

    temp_x = T[]
    temp_y = T[]

    for (xi, yi) in zip(x, y)
        if isnan(xi) || isnan(yi)
            if !isempty(temp_x)
                push!(segments_x, temp_x)
                push!(segments_y, temp_y)
                temp_x = []
                temp_y = []
            end
        else
            push!(temp_x, xi)
            push!(temp_y, yi)
        end
    end

    if !isempty(temp_x)
        push!(segments_x, temp_x)
        push!(segments_y, temp_y)
    end

    return segments_x, segments_y
end

function sign_switch_fractions(vec::AbstractVector)
    # Find sign switches
    sign_switches = findall(i -> sign(vec[i]) != sign(vec[i+1]), 1:length(vec)-1)

    # Calculate interpolation fractions
    fractions = Float64[]
    for i in sign_switches
        Δ = vec[i+1] - vec[i]
        if Δ != 0
            f = -vec[i] / Δ
            push!(fractions, f)
        else
            push!(fractions, NaN)
        end
    end

    # Determine sign of each segment
    segment_signs = []
    prev_idx = 1
    for switch_idx in sign_switches
        push!(segment_signs, sign(vec[prev_idx]))
        prev_idx = switch_idx + 1
    end
    # Add sign of last segment
    push!(segment_signs, sign(vec[prev_idx]))

    return sign_switches, fractions, segment_signs
end

function find_neighbors_within_distance(all_x, all_y, d)
    # Combine vectors into a 2D points matrix
    points = hcat(all_x, all_y)'
    points = ustrip.(Unitful.unit(d), points)

    # Create a k-d tree for fast neighbor search
    tree = KDTree(points)

    # Find neighbors within threshold distance for all points
    indices_within_d = [Int[] for i in eachindex(all_x)]
    dists_within_d = [typeof(d)[] for i in eachindex(all_x)]
    for i in eachindex(all_x)
        idxs, dists = knn(tree, points[:, i], 6, true)
        append!(indices_within_d[i], idxs[dists.≤ustrip(d)])
        append!(dists_within_d[i], dists[dists.≤ustrip(d)] * Unitful.unit(d))
    end

    return indices_within_d, dists_within_d
end