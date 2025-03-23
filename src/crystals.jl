export NonlinearCrystal, UnidirectionalCrystal, BidirectionalCrystal, default_lambda, default_temp, is_lambda_valid, valid_lambda_range, construct_d_tensor, optical_axis_angle, refraction_data_hi_lo, plot_birefringent_refraction, default_temp

abstract type NonlinearCrystal end

function default_lambda(cr::NonlinearCrystal)
    return isnothing(cr.n_x_principal.lambda_range) ? 633u"nm" : sum(cr.n_x_principal.lambda_range) / 2
end

function default_temp(cr::NonlinearCrystal)
    cr.n_x_principal.temp_ref
end

function is_lambda_valid(lambda::Length, cr::NonlinearCrystal)
    return is_lambda_valid(lambda, cr.n_x_principal) && is_lambda_valid(lambda, cr.n_y_principal) && is_lambda_valid(lambda, cr.n_z_principal)
end

function valid_lambda_range(cr::NonlinearCrystal)
    lambda_ranges = [cr.n_x_principal.lambda_range, cr.n_y_principal.lambda_range, cr.n_z_principal.lambda_range]
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

function construct_d_tensor(pointgroup::String, use_kleinman::Bool=true; components...)
    comps = [c[1] for c in components]
    vals = [c[2] for c in components]
    sg = use_kleinman ? symmetry_groups_kleinman : symmetry_groups

    symmetry = get(sg, pointgroup, nothing)
    symmetry !== nothing || error("Point group '$pointgroup' not defined.")

    # Iterate through each symmetry group to ensure consistency
    d_res = zeros(3, 6) * u"pm/V"
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
            d_res[tensor_indices(gc)...] = vals[c_idx] * gs
        end
    end

    return d_res
end

function calc_d_eff(
    cr::NonlinearCrystal, 
    E_dir_r1::AbstractVector{<:Number}, 
    E_dir_r2::AbstractVector{<:Number}, 
    E_dir_b::AbstractVector{<:Number}
)
    P_dir_b = [
        E_dir_r1[1] * E_dir_r2[1],                             # xx
        E_dir_r1[2] * E_dir_r2[2],                             # yy
        E_dir_r1[3] * E_dir_r2[3],                             # zz
        E_dir_r1[2] * E_dir_r2[3] + E_dir_r1[3] * E_dir_r2[2], # yz + zy
        E_dir_r1[1] * E_dir_r2[3] + E_dir_r1[3] * E_dir_r2[1], # xz + zx
        E_dir_r1[1] * E_dir_r2[2] + E_dir_r1[2] * E_dir_r2[1]  # xy + yx
    ]

    return dot(E_dir_b, cr.d * P_dir_b)
end


## Unidirectional crystal

struct UnidirectionalCrystal{TM,TE,TO,TD} <: NonlinearCrystal
    metadata::TM
    n_xy_principal::TE
    n_z_principal::TO
    d::SMatrix{3,6,TD}

    function UnidirectionalCrystal(
        metadata::Dict,
        n_o_principal::RefractiveIndex,
        n_e_principal::RefractiveIndex,
        d::AbstractMatrix,
    )
        return new{typeof(metadata),typeof(n_o_principal),typeof(n_e_principal),eltype(d)}(
            metadata, n_o_principal, n_e_principal, SMatrix{3,6,eltype(d)}(d)
        )
    end
end

function Base.getproperty(cr::UnidirectionalCrystal, sym::Symbol)
    if sym === :n_x_principal || sym === :n_y_principal || sym === :n_o_principal
        return cr.n_xy_principal
    elseif sym === :n_e_principal
        return cr.n_z_principal
    else # Fallback to real fields
        return getfield(cr, sym)
    end
end

## Bidirectional crystal

struct BidirectionalCrystal{TM,TX,TY,TZ,TD} <: NonlinearCrystal
    metadata::TM
    n_x_principal::TX
    n_y_principal::TY
    n_z_principal::TZ
    d::SMatrix{3,6,TD}

    function BidirectionalCrystal(
        metadata::Dict,
        n_x_principal::RefractiveIndex,
        n_y_principal::RefractiveIndex,
        n_z_principal::RefractiveIndex,
        d::AbstractMatrix,
    )
        return new{typeof(metadata),typeof(n_x_principal),typeof(n_y_principal),typeof(n_z_principal),eltype(d)}(
            metadata, n_x_principal, n_y_principal, n_z_principal, SMatrix{3,6,eltype(d)}(d)
        )
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
    o_or_e_r1_r2_b::AbstractVector{Symbol}, 
    cr::UnidirectionalCrystal,
    lambda_r1_r2_b::AbstractVector{<:Unitful.Length};
    temp::Temperature=default_temp(cr),
)
    return [o_e_to_hi_lo(o_or_e_r1_r2_b[i], cr, lambda_r1_r2_b[i]; temp) for i in eachindex(o_or_e_r1_r2_b)]
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
    hi_or_lo_r1_r2_b::AbstractVector{Symbol}, 
    cr::UnidirectionalCrystal,
    lambda_r1_r2_b::AbstractVector{<:Unitful.Length};
    temp::Temperature=default_temp(cr),
)
    return [hi_lo_to_o_e(hi_or_lo_r1_r2_b[i], cr, lambda_r1_r2_b[i]; temp) for i in eachindex(hi_or_lo_r1_r2_b)]
end

function polarization_r1_r2_b_to_hi_lo(
    pol_r1_r2_b::AbstractVector{Symbol}, 
    cr::NonlinearCrystal,
    lambda_r1_r2_b::AbstractVector{<:Unitful.Length};
    temp::Temperature=default_temp(cr),
)
    if all([p in [:hi, :lo] for p in pol_r1_r2_b])
        return pol_r1_r2_b
    elseif isa(cr, UnidirectionalCrystal) 
        if all([p in [:o, :e] for p in pol_r1_r2_b])
            return o_e_to_hi_lo(pol_r1_r2_b, cr, lambda_r1_r2_b; temp)
        else
            @error "$(pol_r1_r2_b) is unvalid polarization data."
        end
    else 
        @error "$(pol_r1_r2_b) is unvalid polarization data for a crystal of type $(typeof(cr))."
    end
end

function optical_axis_angle(
    cr::BidirectionalCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr);
)
    if isnothing(lambda)
        lambda = isnothing(cr.n_x_principal.lambda_range) ? 633u"nm" : sum(cr.n_x_principal.lambda_range) / 2
    end

    nx = refractive_index(cr.n_x_principal, lambda, temp)
    ny = refractive_index(cr.n_y_principal, lambda, temp)
    nz = refractive_index(cr.n_z_principal, lambda, temp)

    # TODO: Clean up using matrix notation
    if nx < nz
        Vz = asin((nz * sqrt(ny^2 - nx^2)) / (ny * sqrt(nz^2 - nx^2)))
    else
        Vz = acos((nx * sqrt(ny^2 - nz^2)) / (ny * sqrt(nx^2 - nz^2)))
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

function refraction_data_hi_lo(
    θ::Angle,
    ϕ::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr);
    n_hi_lo_only::Bool=false,
)
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"

    k_dir, ε_tensor, n_xyz = calc_k_dir_ε_tensor_n_xyz(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    n_hi_lo_only && return n_hi_lo

    E_dir_hi_lo, S_dir_hi_lo = calc_E_dir_S_dir(k_dir, ε_tensor, D_dir_hi_lo)

    walkoff_angle_hi_lo = (acos(clamp(dot(S_dir_hi_lo[1], k_dir), -1, 1)), acos(clamp(dot(S_dir_hi_lo[2], k_dir), -1, 1))) .|> u"rad"

    # Calculate derivative based data
    group_index_hi_lo = calc_group_index_hi_lo(θ, ϕ, cr, lambda, temp)
    β0_hi_lo = calc_β0_hi_lo(θ, ϕ, cr, lambda, temp)
    β1_hi_lo = calc_β1_hi_lo(θ, ϕ, cr, lambda, temp)
    β2_hi_lo = calc_β2_hi_lo(θ, ϕ, cr, lambda, temp)
    β3_hi_lo = calc_β3_hi_lo(θ, ϕ, cr, lambda, temp)

    return n_hi_lo, D_dir_hi_lo, E_dir_hi_lo, S_dir_hi_lo, walkoff_angle_hi_lo, group_index_hi_lo, β0_hi_lo, β1_hi_lo, β2_hi_lo, β3_hi_lo
end

function calc_k_dir_ε_tensor_n_xyz(
    θ::Angle,
    ϕ::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr),
)
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"

    nx = refractive_index(cr.n_x_principal, lambda, temp)
    ny = refractive_index(cr.n_y_principal, lambda, temp)
    nz = refractive_index(cr.n_z_principal, lambda, temp)

    k_dir = angles_to_vector(θ, ϕ)
    ε_tensor = Diagonal(@SVector [1 / nx^2, 1 / ny^2, 1 / nz^2])

    return k_dir, ε_tensor, (nx, ny, nz)
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

    # Calculate Poynting vectors and walkoff angles
    S_dir_hi = normalize(cross(E_dir_hi_lo[1], H_dir_hi_lo[1]))
    S_dir_hi = sign(dot(S_dir_hi, k_dir)) * S_dir_hi # Let S and k point have the same orientation
    S_dir_lo = normalize(cross(E_dir_hi_lo[2], H_dir_hi_lo[2]))
    S_dir_lo = sign(dot(S_dir_lo, k_dir)) * S_dir_lo # Let S and k point have the same orientation
    S_dir_hi_lo = (S_dir_hi, S_dir_lo)

    return E_dir_hi_lo, S_dir_hi_lo
end

# β0 = k
function calc_β0_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    n_hi_lo = refraction_data_hi_lo(θ, ϕ, cr, lambda, temp; n_hi_lo_only=true)
    return Tuple([(2π / lambda * n) |> u"m^-1" for n in n_hi_lo])
end
calc_β0_hi_lo(θ::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β0_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β1 = ∂k/∂ω
function calc_β1_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
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
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
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
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
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
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
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
    k_dir, ε_tensor, (nx, ny, nz) = calc_k_dir_ε_tensor_n_xyz(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo, E_dir_hi_lo, S_dir_hi_lo, walkoff_angle_hi_lo = refraction_data_hi_lo(θ, ϕ, cr, lambda, temp)

    f = Figure()
    ax = Axis3(f[1, 1], azimuth=0.1π, elevation=0.05π, aspect=:data, viewmode=:fit)

    draw_ellipsoid!(ax, [0, 0, 0], [1 / nx^2, 1 / ny^2, 1 / nz^2]; alpha=0.2)

    lines!(ax, [(0, 0, 0), (S_dir_hi_lo[1]...,)], color=1, colormap=:Paired_12, colorrange=(1, 12), label="S_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (S_dir_hi_lo[2]...,)], color=2, colormap=:Paired_12, colorrange=(1, 12), label="S_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (E_dir_hi_lo[1]...,) .* 0.2], color=3, colormap=:Paired_12, colorrange=(1, 12), label="E_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (E_dir_hi_lo[2]...,) .* 0.2], color=4, colormap=:Paired_12, colorrange=(1, 12), label="E_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (D_dir_hi_lo[1]...,) .* 0.3], color=5, colormap=:Paired_12, colorrange=(1, 12), label="D_high", linewidth=5, fxaa=true)
    lines!(ax, [(0, 0, 0), (D_dir_hi_lo[2]...,) .* 0.3], color=6, colormap=:Paired_12, colorrange=(1, 12), label="D_low", linewidth=5, fxaa=true)

    lines!(ax, [(0, 0, 0), (k_dir...,)], color=7, colormap=:Paired_12, colorrange=(1, 12), label="k", linewidth=5, fxaa=true)

    Legend(f[1, 2], ax)

    return f
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

    plot_refractiveindex!(cr.n_x_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_x", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_y_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_y", colormap=Reverse(:Greens))
    plot_refractiveindex!(cr.n_z_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_z", colormap=Reverse(:Blues))
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

    plot_refractiveindex!(cr.n_xy_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_x/n_y (e)", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_z_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_z (o)", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end