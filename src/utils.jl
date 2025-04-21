export all_nonlinearcrystals, all_unidirectionalcrystals, all_bidirectionalcrystals

"""
    subtypes_in_module(mod::Module, T; recurse=false) -> Vector{T}

Returns a list of all objects in module `mod` that are subtypes of type `T`.

If `recurse=true`, submodules will also be searched recursively.
Useful for dynamically collecting types defined under a module namespace.
"""
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

"""
    all_nonlinearcrystals()

Return all available subtypes of `NonlinearCrystal` defined in the `NonlinearCrystals` module.
"""
function all_nonlinearcrystals()
    return subtypes_in_module(NonlinearCrystals, NonlinearCrystal)
end

"""
    all_unidirectionalcrystals()

Return all available subtypes of `UnidirectionalCrystal` defined in the `NonlinearCrystals` module.
"""
function all_unidirectionalcrystals()
    return subtypes_in_module(NonlinearCrystals, UnidirectionalCrystal)
end

"""
    all_bidirectionalcrystals()

Return all available subtypes of `BidirectionalCrystal` defined in the `NonlinearCrystals` module.
"""
function all_bidirectionalcrystals()
    return subtypes_in_module(NonlinearCrystals, BidirectionalCrystal)
end

"""
    auto_fmt(x; digits=3, sci_thresh=1e4, inf_thresh=1e10) -> String

Formats a numeric value `x` into a compact string:
- Uses scientific notation for large/small values
- Displays `"Inf"` or `"NaN"` for extreme cases
- Controls decimal digits via `digits`
"""
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

"""
    vector_to_angles(v::AbstractVector{<:Number}) -> (θ::Angle, ϕ::Angle)

Converts a 3D unit vector `v` to spherical coordinates:
- θ (polar angle): 0 to π
- ϕ (azimuthal angle): normalized to [0, 2π)
"""
function vector_to_angles(v::AbstractVector{<:Number})
    x, y, z = v
    r = sqrt(x^2 + y^2 + z^2)
    θ = Float64(acos(z / r)) |> u"rad" # Polar angle (0 to π)
    ϕ = atan(y, x)  # Azimuthal angle (−π to π)
    ϕ = mod(ϕ, 2π) |> u"rad" # Normalize φ to [0, 2π)
    return θ, ϕ
end

"""
    angles_to_vector(θ::Angle, ϕ::Angle) -> SVector{3,Float64}

Returns a unit direction vector corresponding to spherical angles (θ, ϕ).
"""
function angles_to_vector(theta::Angle, phi::Angle)
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return @SVector [x, y, z]
end

"""
    find_principal_planes(θ::Angle, ϕ::Angle; angle_tol=0.1u"°") -> (Symbol, Symbol)

Returns up to two principal planes (`:XY`, `:XZ`, `:YZ`) to which the vector (θ, ϕ) is approximately aligned.
Principal planes are defined as those orthogonal to the X, Y, or Z axes, within the given angular tolerance `angle_tol`.
"""
function find_principal_planes(theta::Angle, phi::Angle; angle_tol::Angle = 0.1u"°")
    v = angles_to_vector(theta, phi)

    plane_syms = [:YZ, :XZ, :XY]
    plane_normals = [
        (@SVector [1.0, 0.0, 0.0]), 
        (@SVector [0.0, 1.0, 0.0]), 
        (@SVector [0.0, 0.0, 1.0]),
    ]

    p1 = nothing
    p2 = nothing
    count = 0
    for (sym, normal) in zip(plane_syms, plane_normals)
        angle = acos(abs(clamp(dot(v, normal), -1, 1))) |> u"rad"
        if abs(angle - π/2) ≤ angle_tol
            count += 1
            if count == 1
                p1 = sym
            elseif count == 2
                p2 = sym
                # More than two pricipal plane at once is not possible so we stop here
                break  
            end
        end
    end

    return (p1, p2)
end

"""
    axes_to_θ_ϕ(axes::Union{Symbol, Vector{Symbol}}) -> Vector{Tuple{Angle, Angle}}

Maps symbolic crystal axes (`:X`, `:Y`, `:Z`) to corresponding spherical direction angles (θ, ϕ).
"""
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

"""
    bool_permutations(val1, val2, n) -> Vector{NTuple{n, Any}}

Returns all `2^n` permutations of `val1` and `val2` in an `n`-element tuple.
Useful for enumerating polarization combinations (e.g., `:hi` / `:lo`).
"""
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

function eigenvals_2d(A::AbstractMatrix)
    (a, b, c, d) = A'
    @assert ((a + d) / 2)^2 - (a * d - b * c) >= 0.0 "This function cannot handle imaginary eigenvalues."
    discr = sqrt(((a + d) / 2)^2 - (a * d - b * c)) 
    eigvals = (a + d) / 2 .+ (@SVector [-discr, discr])
    return eigvals
end

"""
    eigen_2d(A::AbstractMatrix) -> (eigenvalues::Vector, eigenvectors::SMatrix{2,2})

Returns sorted real eigenvalues and their corresponding orthonormal eigenvectors for a 2×2 matrix `A`.
In comparison with the eigen implementation in Julia's base, this minimal implementation is compatible with ForwardDiff.jl.
"""
function eigen_2d(A::AbstractMatrix)
    eigvals = eigenvals_2d(A')
    
    (a, b, c, d) = A'
    if b == 0 && c == 0 # Catch degenerate case of diagonal matrices
        if a >= d
            return [a, d], @SMatrix [1.0 0.0; 0.0 1.0]
        else
            return [d, a], @SMatrix [0.0 1.0; 1.0 0.0]
        end
    else
        eigvecs_orth1 = [(@SVector [b, eigvals[1] - a]), (@SVector [eigvals[1] - d, c])]
        idx_ev1 = findmax(norm.(eigvecs_orth1))[2]
        eigvecs_orth2 = [(@SVector [b, eigvals[2] - a]), (@SVector [eigvals[2] - d, c])]
        idx_ev2 = findmax(norm.(eigvecs_orth2))[2]
        eigvecs_orth = [normalize(eigvecs_orth1[idx_ev1]) normalize(eigvecs_orth2[idx_ev2])]
        return eigvals, eigvecs_orth
    end
end

"""
    shift_lambda_with_freq(lambda::Length, Δω::Frequency) -> Length

Applies a frequency shift `Δω` to the wavelength `λ`, using the relation:
    ω = 2π·c / λ, and then λ' = 2π·c / (ω + Δω)
Returns the shifted wavelength.
"""
function shift_lambda_with_freq(lambda::Length, delta_omega::Frequency)
    omega = c_0 ./ (2π * lambda)
    lambda_shifted = c_0 / (2π * (omega + delta_omega))
    return lambda_shifted
end

"""
    split_on_nan(x::Vector, y::Vector) -> (Vector{Vector}, Vector{Vector})

Splits two paired vectors into subsegments at `NaN` positions. 
Each returned segment pair contains only valid (non-NaN) values.
"""
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

"""
    sign_switch_fractions(vec::Vector) -> (switch_indices, fractions, segment_signs)

Given a vector of values, detects all zero crossings (i.e., sign changes), and returns:
- indices of sign flips
- linear interpolation fractions where zero-crossings occur
- sign of each segment between crossings
"""
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

"""
    find_neighbors_within_distance(all_x, all_y, d) -> (indices_within_d, dists_within_d)

Given two vectors `all_x` and `all_y` representing 2D point coordinates, and a scalar distance `d` (with units),
this function identifies all neighboring points within distance `d` of each point using a k-d tree.
"""
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