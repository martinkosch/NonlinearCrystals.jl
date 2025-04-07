export NonlinearCrystal, UnidirectionalCrystal, BidirectionalCrystal, default_lambda, default_temp, is_lambda_valid, valid_lambda_range, optical_axis_angle, RefractionDataHiLo, RefractionData, plot_birefringent_refraction, default_temp

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

## Unidirectional crystal

struct UnidirectionalCrystal{TM<:Dict,TE<:RefractiveIndex,TO<:RefractiveIndex} <: NonlinearCrystal
    metadata::TM
    n_XY_principal::TE
    n_Z_principal::TO
    d_XYZ_ref_full::SArray{Tuple{3,3,3},typeof(1.0u"m/V")}
    miller_delta::Union{Nothing,SArray{Tuple{3,3,3},typeof(1.0u"m/V")}}

    function UnidirectionalCrystal(
        metadata::Dict,
        n_o_principal::RefractiveIndex,
        n_e_principal::RefractiveIndex,
        d_XYZ_ref_full::AbstractArray{<:Number,3};
        miller_delta::Union{Nothing,AbstractArray{<:Number,3}}=nothing,
    )
        @assert haskey(metadata, :point_group) && haskey(POINT_GROUP_MAP, metadata[:point_group]) "Point group unknown or not specified in metadata. Possible values are:\n$(keys(NonlinearCrystals.POINT_GROUP_MAP))"
        return new{typeof(metadata),typeof(n_o_principal),typeof(n_e_principal)}(
            metadata,
            n_o_principal,
            n_e_principal,
            d_XYZ_ref_full,
            miller_delta,
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
struct BidirectionalCrystal{TM<:Dict,TX<:RefractiveIndex,TY<:RefractiveIndex,TZ<:RefractiveIndex} <: NonlinearCrystal
    metadata::TM
    n_X_principal::TX
    n_Y_principal::TY
    n_Z_principal::TZ
    d_XYZ_ref_full::SArray{Tuple{3,3,3},typeof(1.0u"m/V")}
    miller_delta::Union{Nothing,SArray{Tuple{3,3,3},typeof(1.0u"m/V")}}

    function BidirectionalCrystal(
        metadata::Dict,
        n_X_principal::RefractiveIndex,
        n_Y_principal::RefractiveIndex,
        n_Z_principal::RefractiveIndex,
        d_XYZ_ref_full::AbstractArray{<:Number,3};
        miller_delta::Union{Nothing,AbstractArray{<:Number,3}}=nothing,
        warn_n_Z_smaller_n_X::Bool=true,
    )
        @assert haskey(metadata, :point_group) && haskey(POINT_GROUP_MAP, metadata[:point_group]) "Point group unknown or not specified in metadata. Possible values are:\n$(keys(NonlinearCrystals.POINT_GROUP_MAP))"
        warnstring = "n_Z_principal ≥ n_Y_principal ≥ n_X_principal should be used in this package for biaxial crystals$(haskey(metadata, :description) ? " but this is not fulfilled for crystal '$(metadata[:description])'." : ".")"
        warn_n_Z_smaller_n_X && (n_Z_principal() < n_Y_principal() < n_X_principal()) && (@warn warnstring)
        return new{typeof(metadata),typeof(n_X_principal),typeof(n_Y_principal),typeof(n_Z_principal)}(
            metadata,
            n_X_principal,
            n_Y_principal,
            n_Z_principal,
            d_XYZ_ref_full,
            miller_delta,
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

function assign_o_or_e(
    principal_plane::Symbol,
    E_dir::AbstractVector{<:Number};
    angle_tol_ud::Angle=0.2u"°",
)
    # For UnidirectionalCrystal, all waves with an E field perpendicular to the optical axis z are ordinary (:o)
    if principal_plane === :UD
        z_vec = @SVector [0.0, 0.0, 1.0]
        angle = acos(clamp(abs(dot(z_vec, E_dir)), -1, 1)) |> u"rad"
        is_o = abs(angle - 90u"°") < angle_tol_ud # Use tolerance due to singularity when propagating close to the optical axis
    else
        # For BidirectionalCrystal, ordinary and extraordinary can only be assigned when propagating within the principal planes
        if principal_plane === :YZ
            plane_normal = @SVector [1.0, 0.0, 0.0]
        elseif principal_plane === :XZ
            plane_normal = @SVector [0.0, 1.0, 0.0]
        elseif principal_plane === :XY
            plane_normal = @SVector [0.0, 0.0, 1.0]
        else
            @error "Principal plane '$(principal_plane)' unknown."
        end

        angle = acos(clamp(abs(dot(plane_normal, E_dir)), -1, 1)) |> u"rad"
        dir_error_str = "E_dir should be either within the plane or parallel to it for stable birefringent modes but it seems to be inbetween"
        @assert (abs(angle - 90u"°") < 5u"°") || (abs(angle) < 5u"°") dir_error_str
        is_o = abs(angle) < 45u"°"
    end


    return is_o ? (:o) : (:e)
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

struct RefractionTypeHiLo
    principal_plane::Union{Symbol}
    o_or_e_hi_lo::NTuple{2,Symbol}
end

function RefractionTypeHiLo(
    principal_plane::Symbol,
    E_dir_hi_lo::NTuple{2,<:AbstractVector{<:Number}};
)
    o_or_e_hi_lo = Tuple(assign_o_or_e(principal_plane, E_dir) for E_dir in E_dir_hi_lo)

    return RefractionTypeHiLo(principal_plane, o_or_e_hi_lo)
end

struct RefractionDataHiLo{CT<:NonlinearCrystal}
    theta::typeof(1.0u"rad")
    phi::typeof(1.0u"rad")
    cr::CT
    lambda::typeof(1.0u"m")
    temp::typeof(1.0u"K")
    refr_type_hi_lo::NTuple{2,Union{Nothing,RefractionTypeHiLo}}
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
            "Refractive index type:",
            "hi ($(rd.refr_type_hi_lo[1].o_or_e_hi_lo[1]))",
            "lo ($(rd.refr_type_hi_lo[1].o_or_e_hi_lo[2]))"
        )
    else
        @printf(io, "%-25s %-25s %-25s\n", "Refractive index type:", "hi", "lo")
        for t in rd.refr_type_hi_lo
            isnothing(t) && continue
            pols = ["$(t.o_or_e_hi_lo[i])" for i in 1:2]
            @printf(io, "%-25s %-25s %-25s\n", "$(t.principal_plane) plane polarization:", pols...)
        end
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
    angle_tol::Angle=0.1u"°",
)
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(θ, ϕ, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)

    E_dir_hi_lo, S_dir_hi_lo = calc_E_dir_S_dir(k_dir, ε_tensor, D_dir_hi_lo)

    if isa(cr, UnidirectionalCrystal)
        # o and e are valid classifiers for UnidirectionalCrystal even if propagating outside of the pricipal planes
        refr_type_hi_lo = (RefractionTypeHiLo(:UD, E_dir_hi_lo), nothing)
    else
        principal_planes = find_principal_planes(θ, ϕ; angle_tol)
        refr_type_hi_lo = Tuple((isnothing(p) ? nothing : RefractionTypeHiLo(p, E_dir_hi_lo) for p in principal_planes))
    end

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
        refr_type_hi_lo,
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
    smallest_axis = findmin(abs.(k_dir))[2]
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

struct RefractionType
    principal_plane::Union{Symbol}
    o_or_e::Symbol
end

function RefractionType(hi_or_lo::Symbol, rt::RefractionTypeHiLo)
    return RefractionType(rt.principal_plane, rt.o_or_e_hi_lo[hi_or_lo === :hi ? 1 : 2])
end

struct RefractionData{CT<:NonlinearCrystal}
    hi_or_lo::Symbol
    theta::typeof(1.0u"rad")
    phi::typeof(1.0u"rad")
    cr::CT
    lambda::typeof(1.0u"m")
    temp::typeof(1.0u"K")
    refr_type::NTuple{2,Union{Nothing,RefractionType}}
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
        @printf(io, "%-25s %-25s\n", "Refractive index type:", string(rd.hi_or_lo) * " (" * string(rd.refr_type[1].o_or_e) * ")")
    else
        @printf(io, "%-25s %-25s\n", "Refractive index type:", rd.hi_or_lo === :hi ? "hi" : "lo")
        for t in rd.refr_type
            isnothing(t) && continue
            @printf(io, "%-25s %-25s\n", "$(t.principal_plane) plane polarization:", "$(t.o_or_e)")
        end
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
    digits = 3

    @printf(io, "%-25s %s\n", "Crystal:", rd.cr.metadata[:description])

    λ = ustrip.(u"nm", round.(u"nm", rd.lambda; digits))
    @printf(io, "%-25s %s\n", "Wavelength (nm):", λ)

    @printf(io, "%-25s θ: %3.2f°, ϕ: %3.2f°\n", "k angles:",
        ustrip(u"°", rd.theta), ustrip(u"°", rd.phi))

    @printf(io, "%-25s %-25s\n", "k direction:",
        vec_str(angles_to_vector(rd.theta, rd.phi))
    )

    @printf(io, "%-25s %3.2f K (%3.2f °C)\n", "Temperature:",
        ustrip(u"K", rd.temp), ustrip(u"°C", rd.temp))
end

function RefractionData(hi_or_lo::Symbol, rd::RefractionDataHiLo{CT}) where {CT}
    fields = map(fieldnames(RefractionDataHiLo)) do f
        f in (:theta, :phi, :cr, :lambda, :temp) && return getfield(rd, f)
        (f === :refr_type_hi_lo) && return Tuple(
            isnothing(p) ? nothing : RefractionType(hi_or_lo, p) for p in getfield(rd, f)
        )
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

    plot_refractiveindex!(cr.n_XY_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_X/n_Y (o)", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_Z_principal; n_sample_pts, temp, lambda_min, lambda_max, label="n_Z (e)", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end