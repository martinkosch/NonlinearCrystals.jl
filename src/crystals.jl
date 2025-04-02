export NonlinearCrystal, UnidirectionalCrystal, BidirectionalCrystal, default_lambda, default_temp, is_lambda_valid, valid_lambda_range, calc_d_XYZ_full, optical_axis_angle, RefractionDataHiLo, RefractionData, plot_birefringent_refraction, default_temp

abstract type NonlinearCrystal end

Base.broadcastable(cr::NonlinearCrystal) = Ref(cr)

function default_lambda(cr::NonlinearCrystal)
    return isnothing(cr.n_X_principal.lambda_range) ? 633u"nm" : sum(cr.n_X_principal.lambda_range) / 2
end

function default_temp(cr::NonlinearCrystal)
    cr.n_X_principal.temp_ref
end

function is_lambda_valid(lambda::Length, cr::NonlinearCrystal; warn_tol::Length=1u"nm")
    return is_lambda_valid(lambda, cr.n_X_principal; warn_tol) && is_lambda_valid(lambda, cr.n_Y_principal; warn_tol) && is_lambda_valid(lambda, cr.n_Z_principal; warn_tol)
end

function valid_lambda_range(cr::NonlinearCrystal)
    lambda_ranges = [cr.n_X_principal.lambda_range, cr.n_Y_principal.lambda_range, cr.n_Z_principal.lambda_range]
    min_valid_lambda = maximum([l[1] for l in lambda_ranges])
    max_valid_lambda = minimum([l[2] for l in lambda_ranges])
    return (min_valid_lambda, max_valid_lambda)
end

function tensor_indices(comp::Symbol)
    s = string(comp)
    length(s) == 3 && s[1] == 'd' || error("Invalid tensor component symbol: $comp")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    (i, j)
end

# Expand 3x6 contracted tensor into full 3x3x3
function expand_d_contract(d_contract::AbstractMatrix{<:Number})
    d_full = zeros(eltype(d_contract), 3, 3, 3)
    idx_map = [
        (1, 1), (2, 2), (3, 3),  # J=1,2,3
        (2, 3), (1, 3), (1, 2)   # J=4,5,6 (symmetric indices)
    ]
    for i in 1:3, J in 1:6
        j, k = idx_map[J]
        d_full[i, j, k] += d_contract[i, J]
        if j != k
            d_full[i, k, j] += d_contract[i, J]  # symmetry
        end
    end
    return d_full
end

# Rotate 3x3x3 tensor
function rotate_tensor3(d::AbstractArray{<:Number,3}, R::AbstractMatrix{<:Number})
    @tullio d_rot[i_dash, j_dash, k_dash] := R[i_dash, i] * R[j_dash, j] * R[k_dash, k] * d[i, j, k]
end

# Contract back to 3x6
function contract_d_tensor(d_full::AbstractArray{<:Number,3})
    d_contract = zeros(eltype(d_full), 3, 6)
    idx_map = [
        (1, 1), (2, 2), (3, 3),
        (2, 3), (1, 3), (1, 2)
    ]
    for i in 1:3, J in 1:6
        j, k = idx_map[J]
        if j == k
            d_contract[i, J] = d_full[i, j, k]
        else
            d_contract[i, J] = 0.5 * (d_full[i, j, k] + d_full[i, k, j])
        end
    end
    return d_contract
end

function rot_mat_abc_to_XYZ(axes_assignment_XYZ::NTuple{3,Symbol})
    # Standard basis in XYZ
    basis_vectors = Dict(
        :X => [1.0, 0.0, 0.0],
        :Y => [0.0, 1.0, 0.0],
        :Z => [0.0, 0.0, 1.0]
    )

    # Create a::Vector, b::Vector, c::Vector based on XYZ assignment
    # We want to know how a, b, c are expressed in XYZ
    axis_order = [:a, :b, :c]
    rot_mat = zeros(3, 3)
    for (i, axis_sym) in enumerate(axis_order)
        # Find which XYZ direction this axis is assigned to
        idx_in_XYZ = findfirst(x -> x == axis_sym, axes_assignment_XYZ)
        rot_mat[i, :] .= basis_vectors[[:X, :Y, :Z][idx_in_XYZ]]
    end

    return rot_mat
end

function calc_d_XYZ_full(
    point_group::String,
    rot_mat::AbstractMatrix{<:Number}=I(3),
    use_kleinman::Bool=true;
    components_abc...
)
    comps = [c[1] for c in components_abc]
    vals = [c[2] for c in components_abc]
    sg = use_kleinman ? SYMMETRY_GROUPS_KLEINMAN : SYMMETRY_GROUPS

    symmetry = get(sg, point_group, nothing)
    symmetry !== nothing || error("Point group '$point_group' not defined.")

    # Iterate through each symmetry group to ensure consistency
    # Symmetry goups and given components are assumed to be specified in abc coordinates
    d_abc = zeros(3, 6) * u"pm/V"
    for group in symmetry
        group_comps, group_signs = group
        given_group_comps = intersect(Set(group_comps), Set(comps))

        if length(given_group_comps) > 1
            @error "Exactly one of the following components must be specified: $(group_comps)"
        elseif length(given_group_comps) < 1
            if length(group_comps) == 1
                @error "Component '$(group_comps[1])' must be specified."
            else
                @error "One of these components must be specified: $(group_comps)"
            end
        end

        # Correct signs within the current symmetry group
        given_group_comp = collect(given_group_comps)[1]
        idx = findall(c -> (c === given_group_comp), group_comps)[1]
        group_signs_switched = group_signs[idx] == 1 ? group_signs : -group_signs

        # Set all group components based on the given component with the correct sign relations
        c_idx = findall(c -> (c === given_group_comp), comps)[1]
        for (gc, gs) in zip(group_comps, group_signs_switched)
            d_abc[tensor_indices(gc)...] = vals[c_idx] * gs
        end
    end

    # Rotate d tensor from the coordinates system of the given d tensor components to XYZ coordinates based on the provided assigment
    d_abc_full = expand_d_contract(d_abc)
    d_XYZ_full = rotate_tensor3(d_abc_full, rot_mat)
    return d_XYZ_full
end

function calc_d_eff(
    cr::NonlinearCrystal,
    E_dir_r1::AbstractVector{<:Number},
    E_dir_r2::AbstractVector{<:Number},
    E_dir_b::AbstractVector{<:Number};
    lambda_rrb::Union{AbstractVector{<:Length},Nothing}=nothing,
    temp::Temperature=default_temp(cr),
)
    # If all lambdas are given, use Miller scaling
    # if isnothing(lambda_rrb)
    @tullio d_eff := E_dir_b[i] * cr.d_XYZ_ref_full[i, j, k] * E_dir_r1[j] * E_dir_r2[k]
    # else
    #     # Compute refractive indices along XYZ at given wavelengths and temperature
    #     n_r1 = @SVector [
    #         cr.n_X_principal(lambda_rrb[1], temp),
    #         cr.n_Y_principal(lambda_rrb[1], temp),
    #         cr.n_Z_principal(lambda_rrb[1], temp)
    #     ]
    #     n_r2 = @SVector [
    #         cr.n_X_principal(lambda_rrb[2], temp),
    #         cr.n_Y_principal(lambda_rrb[2], temp),
    #         cr.n_Z_principal(lambda_rrb[2], temp)
    #     ]
    #     n_b = @SVector [
    #         cr.n_X_principal(lambda_rrb[3], temp),
    #         cr.n_Y_principal(lambda_rrb[3], temp),
    #         cr.n_Z_principal(lambda_rrb[3], temp)
    #     ]

    #     # Use Miller-scaled d_XYZ tensor to calculate d_eff
    #     @tullio d_eff := E_dir_b[i] * cr.miller_delta[i, j, k] * n_r1[i] * n_r2[j] * n_b[k] * E_dir_r1[j] * E_dir_r2[k]
    # end
    return d_eff
end

# function compute_miller_delta(
#     d_ref_XYZ::AbstractMatrix{<:Number},
#     lambda_r1::Length,
#     lambda_r2::Length,
#     lambda_b::Length,
#     temp::Temperature,
#     n_X_principal::RefractiveIndex,
#     n_Y_principal::RefractiveIndex,
#     n_Z_principal::RefractiveIndex,
# )
#     @assert size(d_ref_XYZ) == (3, 6)
#     d_ref_XYZ_expanded = expand_d_contract(d_ref_XYZ)

#     n_r1 = @SVector [n_X_principal(lambda_r1, temp), n_Y_principal(lambda_r1, temp), n_Z_principal(lambda_r1, temp)]
#     n_r2 = @SVector [n_X_principal(lambda_r2, temp), n_Y_principal(lambda_r2, temp), n_Z_principal(lambda_r2, temp)]
#     n_b = @SVector [n_X_principal(lambda_b, temp), n_Y_principal(lambda_b, temp), n_Z_principal(lambda_b, temp)]

#     @tullio miller_delta[i, j, k] := d_ref_XYZ_expanded[i, j, k] / (n_r1[i] * n_r2[j] * n_b[k])

#     return miller_delta
# end

## Unidirectional crystal

struct UnidirectionalCrystal{TM,TE,TO,TD} <: NonlinearCrystal
    metadata::TM
    n_XY_principal::TE
    n_Z_principal::TO
    d_XYZ_ref_full::SArray{Tuple{3,3,3},TD}

    function UnidirectionalCrystal(
        metadata::Dict,
        n_o_principal::RefractiveIndex,
        n_e_principal::RefractiveIndex,
        d_XYZ_ref_full::AbstractArray{<:Number,3};
    )
        @assert haskey(metadata, :point_group) && haskey(POINT_GROUP_MAP, metadata[:point_group]) "Point group unknown or not specified in metadata. Possible values are:\n$(keys(NonlinearCrystals.POINT_GROUP_MAP))"
        return new{typeof(metadata),typeof(n_o_principal),typeof(n_e_principal),eltype(d_XYZ_ref_full)}(
            metadata, n_o_principal, n_e_principal, SArray{Tuple{3,3,3},eltype(d_XYZ_ref_full)}(d_XYZ_ref_full)
        )
    end
end

function Base.getproperty(cr::UnidirectionalCrystal, sym::Symbol)
    if sym === :n_X_principal || sym === :n_Y_principal || sym === :n_o_principal
        return cr.n_XY_principal
    elseif sym === :n_e_principal
        return cr.n_Z_principal
    elseif sym === :d_XYZ_ref
        return contract_d_tensor(cr.d_XYZ_ref_full)
    else # Fallback to real fields
        return getfield(cr, sym)
    end
end

## Bidirectional crystal

struct BidirectionalCrystal{TM,TX,TY,TZ,TD} <: NonlinearCrystal
    metadata::TM
    n_X_principal::TX
    n_Y_principal::TY
    n_Z_principal::TZ
    d_XYZ_ref_full::SArray{Tuple{3,3,3},TD}

    function BidirectionalCrystal(
        metadata::Dict,
        n_X_principal::RefractiveIndex,
        n_Y_principal::RefractiveIndex,
        n_Z_principal::RefractiveIndex,
        d_XYZ_ref_full::AbstractArray{<:Number,3};
        warn_n_Z_smaller_n_X::Bool=true,
    )
        @assert haskey(metadata, :point_group) && haskey(POINT_GROUP_MAP, metadata[:point_group]) "Point group unknown or not specified in metadata. Possible values are:\n$(keys(NonlinearCrystals.POINT_GROUP_MAP))"
        warnstring = "n_Z_principal ≥ n_Y_principal ≥ n_X_principal should be used in this package for biaxial crystals$(haskey(metadata, :description) ? " but this is not fulfilled for crystal '$(metadata[:description])'." : ".")"
        warn_n_Z_smaller_n_X && (n_Z_principal() < n_Y_principal() < n_X_principal()) && (@warn warnstring)
        return new{typeof(metadata),typeof(n_X_principal),typeof(n_Y_principal),typeof(n_Z_principal),eltype(d_XYZ_ref_full)}(
            metadata, n_X_principal, n_Y_principal, n_Z_principal, SArray{Tuple{3,3,3},eltype(d_XYZ_ref_full)}(d_XYZ_ref_full)
        )
    end
end

function Base.getproperty(cr::BidirectionalCrystal, sym::Symbol)
    if sym === :d_XYZ_ref
        return contract_d_tensor(cr.d_XYZ_ref_full)
    else # Fallback to real fields
        return getfield(cr, sym)
    end
end

function o_e_to_hi_lo(
    o_or_e::Symbol,
    cr::UnidirectionalCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
)
    n_o_principal = refractive_index(cr.n_o_principal, lambda, temp)
    n_e_principal = refractive_index(cr.n_e_principal, lambda, temp)

    o_or_e === :o && return n_o_principal > n_e_principal ? (:hi) : (:lo)
    o_or_e === :e && return n_e_principal > n_o_principal ? (:hi) : (:lo)
    error("Polarization must be :o or :e. Currently it is $(o_or_e).")
end

function o_e_to_hi_lo(
    o_or_e_rrb::NTuple{3,Symbol},
    cr::UnidirectionalCrystal,
    lambda_rrb::NTuple{3,Length};
    temp::Temperature=default_temp(cr),
)
    return Tuple(o_e_to_hi_lo(o_or_e_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(o_or_e_rrb))
end

function hi_lo_to_o_e(
    hi_or_lo::Symbol,
    cr::UnidirectionalCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
)
    n_o_principal = refractive_index(cr.n_o_principal, lambda, temp)
    n_e_principal = refractive_index(cr.n_e_principal, lambda, temp)

    hi_or_lo === :hi && return n_o_principal > n_e_principal ? :o : :e
    hi_or_lo === :lo && return n_o_principal > n_e_principal ? :e : :o
    error("Polarization must be :hi or :lo. Currently it is '$(hi_or_lo)'.")
end

function hi_lo_to_o_e(
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::UnidirectionalCrystal,
    lambda_rrb::NTuple{3,Length};
    temp::Temperature=default_temp(cr),
)
    return Tuple(hi_lo_to_o_e(hi_or_lo_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(hi_or_lo_rrb))
end

function polarization_rrb_to_hi_lo(
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal,
    lambda_rrb::NTuple{3,Length};
    temp::Temperature=default_temp(cr),
)
    if all([p in [:hi, :lo] for p in pol_rrb])
        return pol_rrb
    elseif isa(cr, UnidirectionalCrystal)
        if all([p in [:o, :e] for p in pol_rrb])
            return o_e_to_hi_lo(pol_rrb, cr, lambda_rrb; temp)
        else
            @error "$(pol_rrb) is unvalid polarization data."
        end
    else
        @error "$(pol_rrb) is unvalid polarization data for a crystal of type $(typeof(cr))."
    end
end

function optical_axis_angle(
    cr::BidirectionalCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr);
)
    if isnothing(lambda)
        lambda = isnothing(cr.n_X_principal.lambda_range) ? 633u"nm" : sum(cr.n_X_principal.lambda_range) / 2
    end

    n_X = refractive_index(cr.n_X_principal, lambda, temp)
    n_Y = refractive_index(cr.n_Y_principal, lambda, temp)
    n_Z = refractive_index(cr.n_Z_principal, lambda, temp)

    # TODO: Clean up using matrix notation
    if n_X < n_Z
        Vz = asin((n_Z * sqrt(n_Y^2 - n_X^2)) / (n_Y * sqrt(n_Z^2 - n_X^2)))
    else
        Vz = acos((n_X * sqrt(n_Y^2 - n_Z^2)) / (n_Y * sqrt(n_X^2 - n_Z^2)))
    end

    return Vz |> u"rad"
end

function optical_axis_angle(
    cr::UnidirectionalCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr),
)
    return 0.0u"rad"
end

struct RefractionDataHiLo{CT<:NonlinearCrystal}
    theta::typeof(1.0u"rad")
    phi::typeof(1.0u"rad")
    cr::CT
    lambda::typeof(1.0u"m")
    temp::typeof(1.0u"K")
    n_hi_lo::NTuple{2,Float64}
    group_index_hi_lo::NTuple{2,Float64}
    walkoff_angle_hi_lo::NTuple{2,typeof(1.0u"rad")}
    D_dir_hi_lo::NTuple{2,SVector{3,Float64}}
    E_dir_hi_lo::NTuple{2,SVector{3,Float64}}
    S_dir_hi_lo::NTuple{2,SVector{3,Float64}}
    beta0_hi_lo::NTuple{2,typeof(1.0u"m^-1")}
    beta1_hi_lo::NTuple{2,typeof(1.0u"s * m^-1")}
    beta2_hi_lo::NTuple{2,typeof(1.0u"s^2 * m^-1")}
    beta3_hi_lo::NTuple{2,typeof(1.0u"s^3 * m^-1")}
end

function Base.show(io::IO, rd::RefractionDataHiLo)
    digits = 3

    # Print header
    print_refraction_data_header(io, rd)

    println(io, "───────────────────────────────────────────────────────────────────────────")

    # Print content
    if isa(rd.cr, UnidirectionalCrystal)
        @printf(io, "%-25s %-25s %-25s\n",
            "Polarization:",
            "hi ($(hi_lo_to_o_e(:hi, rd.cr, rd.lambda; rd.temp)))",
            "lo ($(hi_lo_to_o_e(:lo, rd.cr, rd.lambda; rd.temp)))"
        )
    else
        @printf(io, "%-25s %-25s %-25s\n", "", "hi", "lo")
    end

    rd_hi = RefractionData(:hi, rd)
    rd_lo = RefractionData(:lo, rd)

    @printf(io, "%-25s %-25s %-25s\n", "Phase velocity / c₀:", auto_fmt(rd_hi.n; digits), auto_fmt(rd_lo.n; digits))
    @printf(io, "%-25s %-25s %-25s\n", "Group velocity / c₀:", auto_fmt(rd_hi.group_index), auto_fmt(rd_lo.group_index))
    @printf(io, "%-25s %-25s %-25s\n", "Walkoff angle (mrad):",
        auto_fmt(ustrip(u"mrad", rd_hi.walkoff_angle); digits),
        auto_fmt(ustrip(u"mrad", rd_lo.walkoff_angle); digits),
    )
    @printf(io, "%-25s %-25s %-25s\n", "S direction:", vec_str(rd_hi.S_dir), vec_str(rd_lo.S_dir))
    @printf(io, "%-24s %-25s %-25s\n", "E direction:", "±" * vec_str(rd_hi.E_dir), "±" * vec_str(rd_lo.E_dir))
    @printf(io, "%-24s %-25s %-25s\n", "D direction:", "±" * vec_str(rd_hi.D_dir), "±" * vec_str(rd_lo.D_dir))
    @printf(io, "%-25s %-25s %-25s\n", "GDD (fs²/mm):",
        auto_fmt(ustrip(u"fs^2/mm", rd_hi.beta2); digits),
        auto_fmt(ustrip(u"fs^2/mm", rd_lo.beta2); digits),
    )
    @printf(io, "%-25s %-25s %-25s\n", "TOD (fs³/mm):",
        auto_fmt(ustrip(u"fs^3/mm", rd_hi.beta3); digits),
        auto_fmt(ustrip(u"fs^3/mm", rd_lo.beta3); digits),
    )
end

function RefractionDataHiLo(
    θ::Angle,
    ϕ::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
)
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)

    E_dir_hi_lo, S_dir_hi_lo = calc_E_dir_S_dir(k_dir, ε_tensor, D_dir_hi_lo)

    walkoff_angle_hi_lo = (s -> (acos(clamp(dot(s, k_dir), -1, 1)))).(S_dir_hi_lo) .|> u"rad"

    # Calculate derivative based data
    group_index_hi_lo = calc_group_index_hi_lo(θ, ϕ, cr, lambda, temp)
    β0_hi_lo = calc_β0_hi_lo(θ, ϕ, cr, lambda, temp)
    β1_hi_lo = calc_β1_hi_lo(θ, ϕ, cr, lambda, temp)
    β2_hi_lo = calc_β2_hi_lo(θ, ϕ, cr, lambda, temp)
    β3_hi_lo = calc_β3_hi_lo(θ, ϕ, cr, lambda, temp)

    rd = RefractionDataHiLo{typeof(cr)}(
        θ,
        ϕ,
        cr,
        lambda,
        temp,
        n_hi_lo,
        group_index_hi_lo,
        walkoff_angle_hi_lo,
        D_dir_hi_lo,
        E_dir_hi_lo,
        S_dir_hi_lo,
        β0_hi_lo,
        β1_hi_lo,
        β2_hi_lo,
        β3_hi_lo,
    )

    return rd
end

# Minmal refractive index calculation, suited for derivatives via ForwardDiff.jl
function calc_n_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length; temp::Temperature)
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    return n_hi_lo
end

function calc_k_dir_ε_tensor_n_XYZ(
    θ::Angle,
    ϕ::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr),
)
    n_X = refractive_index(cr.n_X_principal, lambda, temp)
    n_Y = refractive_index(cr.n_Y_principal, lambda, temp)
    n_Z = refractive_index(cr.n_Z_principal, lambda, temp)

    k_dir = angles_to_vector(θ, ϕ)
    ε_tensor = Diagonal(@SVector [1 / n_X^2, 1 / n_Y^2, 1 / n_Z^2])

    return k_dir, ε_tensor, (n_X, n_Y, n_Z)
end

function calc_n_hi_lo_D_dir_hi_lo(k_dir::AbstractVector{<:Real}, ε_tensor::AbstractMatrix{<:Real})
    # Construct numerically stable coordinate system around k   
    smallest_axis = findmin(k_dir)[2]
    k_orth1 = normalize(cross(k_dir, [i == smallest_axis ? 1.0 : 0.0 for i in 1:3]))
    k_orth2 = normalize(cross(k_dir, k_orth1))
    P_k_orth = [k_orth1 k_orth2]

    # Search for eigenvalues and vectors in the reduced P_k_orth 2D space
    ε_proj = P_k_orth' * ε_tensor * P_k_orth
    eigvals, eigvecs_orth = eigen_2d(ε_proj')
    eigvecs = P_k_orth * eigvecs_orth

    idx = sortperm(eigvals, rev=false)
    n_hi_lo = (1 / sqrt(eigvals[idx[1]]), 1 / sqrt(eigvals[idx[2]]))
    D_dir_hi_lo = (eigvecs[:, idx[1]], eigvecs[:, idx[2]])
    return n_hi_lo, D_dir_hi_lo
end

function calc_E_dir_S_dir(
    k_dir::AbstractVector{<:Real},
    ε_tensor::AbstractMatrix{<:Real},
    D_dir_hi_lo::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}},
)
    # Polarization directions
    E_dir_hi_lo = (normalize(ε_tensor \ D_dir_hi_lo[1]), normalize(ε_tensor \ D_dir_hi_lo[2]))
    H_dir_hi_lo = (cross(k_dir, E_dir_hi_lo[1]), cross(k_dir, E_dir_hi_lo[2]))

    # Calculate Poynting vectors and walkoff k polar angles
    S_dir_hi = normalize(cross(E_dir_hi_lo[1], H_dir_hi_lo[1]))
    S_dir_hi = sign(dot(S_dir_hi, k_dir)) * S_dir_hi # Let S and k point have the same orientation
    S_dir_lo = normalize(cross(E_dir_hi_lo[2], H_dir_hi_lo[2]))
    S_dir_lo = sign(dot(S_dir_lo, k_dir)) * S_dir_lo # Let S and k point have the same orientation
    S_dir_hi_lo = (S_dir_hi, S_dir_lo)

    return E_dir_hi_lo, S_dir_hi_lo
end

# β0 = k
function calc_β0_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    return Tuple([(2π / lambda * n) |> u"m^-1" for n in n_hi_lo])
end
calc_β0_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β0_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β1 = ∂k/∂ω
function calc_β1_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β0_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s/m")
end
calc_β1_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β1_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β2 = ∂²k/∂ω²
function calc_β2_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β1_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^2/m")
end
calc_β2_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β2_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β3 = ∂³k/∂ω³
function calc_β3_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β2_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^3/m")
end
calc_β3_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β3_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

function calc_group_index_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    return Tuple([β1 * c_0 for β1 in calc_β1_hi_lo(θ, ϕ, cr, lambda, temp)])
end
calc_group_index_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_group_index_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

function plot_birefringent_refraction(
    θ,
    ϕ,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr);
)
    k_dir, ε_tensor, (n_X, n_Y, n_Z) = calc_k_dir_ε_tensor_n_XYZ(θ, ϕ, cr, lambda, temp)
    rd = RefractionDataHiLo(θ, ϕ, cr, lambda; temp)

    f = Figure()
    ax = Axis3(f[1, 1], azimuth=0.1π, elevation=0.05π, aspect=:data, viewmode=:fit)

    draw_ellipsoid!(ax, [0, 0, 0], [1 / n_X^2, 1 / n_Y^2, 1 / n_Z^2]; alpha=0.2)

    lines!(ax, [(0, 0, 0), (rd.S_dir_hi_lo[1]...,)], color=1, colormap=:Paired_12, colorrange=(1, 12), label="S_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (rd.S_dir_hi_lo[2]...,)], color=2, colormap=:Paired_12, colorrange=(1, 12), label="S_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (rd.E_dir_hi_lo[1]...,) .* 0.2], color=3, colormap=:Paired_12, colorrange=(1, 12), label="E_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (rd.E_dir_hi_lo[2]...,) .* 0.2], color=4, colormap=:Paired_12, colorrange=(1, 12), label="E_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (rd.D_dir_hi_lo[1]...,) .* 0.3], color=5, colormap=:Paired_12, colorrange=(1, 12), label="D_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (rd.D_dir_hi_lo[2]...,) .* 0.3], color=6, colormap=:Paired_12, colorrange=(1, 12), label="D_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (rd.k_dir...,)], color=7, colormap=:Paired_12, colorrange=(1, 12), label="k", linewidth=5, fxaa=true)

    Legend(f[1, 2], ax)

    return f
end

struct RefractionData{CT<:NonlinearCrystal}
    hi_or_lo::Symbol
    theta::typeof(1.0u"rad")
    phi::typeof(1.0u"rad")
    cr::CT
    lambda::typeof(1.0u"m")
    temp::typeof(1.0u"K")
    n::Float64
    group_index::Float64
    walkoff_angle::typeof(1.0u"rad")
    D_dir::SVector{3,Float64}
    E_dir::SVector{3,Float64}
    S_dir::SVector{3,Float64}
    beta0::typeof(1.0u"m^-1")
    beta1::typeof(1.0u"s * m^-1")
    beta2::typeof(1.0u"s^2 * m^-1")
    beta3::typeof(1.0u"s^3 * m^-1")
end

function Base.show(io::IO, rd::RefractionData)
    digits = 3

    # Print header
    print_refraction_data_header(io, rd)

    # Print content
    if isa(rd.cr, UnidirectionalCrystal)
        if rd.hi_or_lo === :hi
            @printf(io, "%-25s %-25s\n",
                "",
                "hi ($(hi_lo_to_o_e(:hi, rd.cr, rd.lambda; rd.temp)))",
            )
        else
            @printf(io, "%-25s %-25s\n",
                "Polarization:",
                "lo ($(hi_lo_to_o_e(:lo, rd.cr, rd.lambda; rd.temp)))",
            )
        end
    else
        @printf(io, "%-25s %-25s\n", "", rd.hi_or_lo === :hi ? "hi" : "lo")
    end

    @printf(io, "%-25s %-25s\n", "Phase velocity / c₀:", auto_fmt(rd.n; digits))
    @printf(io, "%-25s %-25s\n", "Group velocity / c₀:", auto_fmt(rd.group_index))
    @printf(io, "%-25s %-25s\n", "Walkoff angle (mrad):",
        auto_fmt(ustrip(u"mrad", rd.walkoff_angle); digits),
    )
    @printf(io, "%-25s %-25s\n", "S direction:", vec_str(rd.S_dir))
    @printf(io, "%-24s %-25s\n", "E direction:", "±" * vec_str(rd.E_dir))
    @printf(io, "%-24s %-25s\n", "D direction:", "±" * vec_str(rd.D_dir))
    @printf(io, "%-25s %-25s\n", "GDD (fs²/mm):",
        auto_fmt(ustrip(u"fs^2/mm", rd.beta2); digits),
    )
    @printf(io, "%-25s %-25s\n", "TOD (fs³/mm):",
        auto_fmt(ustrip(u"fs^3/mm", rd.beta3); digits),
    )
end

function print_refraction_data_header(io::IO, rd::Union{RefractionDataHiLo,RefractionData})
    @printf(io, "%-25s %s\n", "Crystal:", rd.cr.metadata[:description])

    @printf(io, "%-25s %s\n", "Wavelength:", uconvert(u"µm", rd.lambda))

    @printf(io, "%-25s θ: %3.1f°, ϕ: %3.1f°\n", "k angles:",
        ustrip(u"°", rd.theta), ustrip(u"°", rd.phi))

    @printf(io, "%-25s %-25s\n", "k direction:",
        vec_str(angles_to_vector(rd.theta, rd.phi))
    )

    @printf(io, "%-25s %3.1f K (%3.1f °C)\n", "Temperature:",
        ustrip(u"K", rd.temp), ustrip(u"°C", rd.temp))
end

function RefractionData(hi_or_lo::Symbol, rd::RefractionDataHiLo{CT}) where {CT}
    fields = map(fieldnames(RefractionDataHiLo)) do f
        f in (:theta, :phi, :cr, :lambda, :temp) && return getfield(rd, f)
        return getfield(rd, f)[hi_or_lo == :hi ? 1 : 2]
    end
    return RefractionData{CT}(hi_or_lo, fields...)
end

function RefractionData(
    hi_or_lo::Symbol,
    θ::Angle,
    ϕ::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
)
    return RefractionData(hi_or_lo, RefractionDataHiLo(θ, ϕ, cr, lambda; temp))
end

function plot_refractiveindex(
    cr::BidirectionalCrystal;
    n_sample_pts=500,
    temp::Union{AbstractVector{<:Unitful.Temperature},Unitful.Temperature}=[default_temp(cr)],
    lambda_min=nothing,
    lambda_max=nothing,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")

    plot_refractiveindex!(cr.n_X_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_X", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_Y_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_Y", colormap=Reverse(:Greens))
    plot_refractiveindex!(cr.n_Z_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_Z", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end

function plot_refractiveindex(
    cr::UnidirectionalCrystal;
    n_sample_pts=500,
    temp::Union{AbstractVector{<:Unitful.Temperature},Unitful.Temperature}=[default_temp(cr)],
    lambda_min=nothing,
    lambda_max=nothing,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")

    plot_refractiveindex!(cr.n_XY_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_X/n_Y (e)", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_Z_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_Z (o)", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end