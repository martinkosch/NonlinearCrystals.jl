export NonlinearCrystal, UnidirectionalCrystal, BidirectionalCrystal, refractive_index_o, refractive_index_e, group_index_o, group_index_e, optical_axis_angle, refraction_data_hi_lo, is_negative_crystal, walkoff_angle, k_dir, plot_birefringent_refraction, phase_match_wavelength, find_all_phase_matches

abstract type NonlinearCrystal end

function default_T(cr::NonlinearCrystal)
    cr.n_x_principal.T_ref
end

function default_λ(cr::NonlinearCrystal)
    return isnothing(cr.n_x_principal.lambda_range) ? 633u"nm" : sum(cr.n_x_principal.lambda_range) / 2
end

calc_k_dir(θ, ϕ) = angles_to_vector(θ, ϕ)

function construct_d_tensor(pointgroup::String; kwargs...)
    d = zeros(3, 6) * u"pm/V"

    if pointgroup == "-43m"
        # Allowed coefficients
        expected_keys = [:d14]
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For '-43m', only the following keyword must be provided: $(expected_keys)"

        d14 = kwargs[:d14]
        d[1, 4] = d14  # d_x_yz
        d[2, 5] = d14  # d_y_xz
        d[3, 6] = d14  # d_z_xy

    elseif pointgroup == "3m"
        expected_keys = [:d31, :d33, :d22, :d15]
        # d15 may be zero or optional, but it is exptected to be explicitly provided
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For '3m', keywords required: $(expected_keys)"

        d31 = kwargs[:d31]
        d33 = kwargs[:d33]
        d22 = kwargs[:d22]
        d15 = kwargs[:d15]
        d[1, 5] = d15       # d_x_xz
        d[1, 6] = -d22      # d_x_xy
        d[2, 1] = -d22      # d_y_xx
        d[2, 2] = d22       # d_y_yy
        d[2, 6] = d15       # d_y_xy
        d[3, 1] = d31       # d_z_xx
        d[3, 2] = d31       # d_z_yy
        d[3, 3] = d33       # d_z_zz

    elseif pointgroup == "mm2"
        expected_keys = [:d15, :d24, :d31, :d32, :d33]
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For 'mm2', keywords required: $(expected_keys)"

        d[1, 5] = kwargs[:d15]   # d_x_xz
        d[2, 4] = kwargs[:d24]   # d_y_yz
        d[3, 1] = kwargs[:d31]   # d_z_xx
        d[3, 2] = kwargs[:d32]   # d_z_yy
        d[3, 3] = kwargs[:d33]   # d_z_zz

    elseif pointgroup == "4mm"
        expected_keys = [:d15, :d31, :d33]
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For '4mm', keywords required: $(expected_keys)"

        d15 = kwargs[:d15]
        d31 = kwargs[:d31]
        d33 = kwargs[:d33]
        d[1, 5] = d15      # d_x_xz
        d[2, 4] = d15      # d_y_yz
        d[3, 1] = d31      # d_z_xx
        d[3, 2] = d31      # d_z_yy
        d[3, 3] = d33      # d_z_zz

    elseif pointgroup == "32"
        expected_keys = [:d11]
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For '32', keyword required: $(expected_keys)"

        d11 = kwargs[:d11]
        d[1, 1] = -d11     # d_x_xx
        d[1, 2] = d11      # d_x_yy
        d[2, 6] = -d11     # d_y_xy

    elseif pointgroup == "6mm"
        expected_keys = [:d11]
        @assert all(sort(collect(keys(kwargs))) .== sort(expected_keys)) "For '6mm', keyword required: $(expected_keys)"

        d11 = kwargs[:d11]
        d[1, 6] = -d11     # d_x_xy
        d[2, 1] = -d11     # d_y_xx
        d[2, 2] = d11      # d_y_yy

    elseif pointgroup == "m-3m"
        @assert isempty(kwargs) "For centrosymmetric 'm-3m', no nonlinear coefficients allowed!"

    else
        error("Unsupported or unknown point group: $pointgroup")
    end

    return d
end

function calc_d_eff(cr::NonlinearCrystal, E_dir_r1::AbstractVector{<:Real}, E_dir_r2::AbstractVector{<:Real}, E_dir_b::AbstractVector{<:Real})
    # Build contracted input combinations:
    P = [
        E_dir_r1[1] * E_dir_r2[1],                             # xx
        E_dir_r1[2] * E_dir_r2[2],                             # yy
        E_dir_r1[3] * E_dir_r2[3],                             # zz
        E_dir_r1[2] * E_dir_r2[3] + E_dir_r1[3] * E_dir_r2[2], # yz + zy
        E_dir_r1[1] * E_dir_r2[3] + E_dir_r1[3] * E_dir_r2[1], # xz + zx
        E_dir_r1[1] * E_dir_r2[2] + E_dir_r1[2] * E_dir_r2[1]  # xy + yx
    ]

    # Contract with blue beam polarization:
    return dot(E_dir_b, cr.d * P)
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


function optical_axis_angle(
    cr::BidirectionalCrystal,
    λ::Unitful.Length=default_λ(cr),
    T::Unitful.Temperature=default_T(cr);
)
    if isnothing(λ)
        λ = isnothing(cr.n_x_principal.lambda_range) ? 633u"nm" : sum(cr.n_x_principal.lambda_range) / 2
    end

    nx = refractive_index(cr.n_x_principal, λ, T)
    ny = refractive_index(cr.n_y_principal, λ, T)
    nz = refractive_index(cr.n_z_principal, λ, T)

    if nx < nz
        Vz = asin((nz * sqrt(ny^2 - nx^2)) / (ny * sqrt(nz^2 - nx^2)))
    else
        Vz = acos((nx * sqrt(ny^2 - nz^2)) / (ny * sqrt(nx^2 - nz^2)))
    end

    return Vz |> u"rad"
end

function optical_axis_angle(
    cr::UnidirectionalCrystal,
    λ::Unitful.Length=default_λ(cr),
    T::Unitful.Temperature=default_T(cr),
)
    return 0.0u"rad"
end

function refraction_data_hi_lo(
    θ,
    ϕ,
    cr::NonlinearCrystal,
    λ::Unitful.Length=default_λ(cr),
    T::Unitful.Temperature=default_T(cr);
    n_hi_lo_only::Bool=false,
)
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"

    k_dir, ε_tensor, n_xyz = calc_k_dir_ε_tensor_n_xyz(θ, ϕ, cr, λ, T)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    n_hi_lo_only && return n_hi_lo

    E_dir_hi_lo, S_dir_hi_lo = calc_E_dir_S_dir(k_dir, ε_tensor, D_dir_hi_lo)

    walkoff_angle_hi_lo = (acos(clamp(dot(S_dir_hi_lo[1], k_dir), -1, 1)), acos(clamp(dot(S_dir_hi_lo[2], k_dir), -1, 1))) .|> u"rad"

    # Assign ordinary or extraordinary direction labels to the calculated hi and lo direction in case of unidirectional crystal
    o_or_e_hi_lo = assign_e_o(E_dir_hi_lo)

    #Calculate derivative based data
    group_index_hi_lo = calc_group_index_hi_lo(θ, ϕ, cr, λ, T)
    β0_hi_lo = calc_β0_hi_lo(θ, ϕ, cr, λ, T)
    β1_hi_lo = calc_β1_hi_lo(θ, ϕ, cr, λ, T)
    β2_hi_lo = calc_β2_hi_lo(θ, ϕ, cr, λ, T)
    β3_hi_lo = calc_β3_hi_lo(θ, ϕ, cr, λ, T)

    return n_hi_lo, D_dir_hi_lo, E_dir_hi_lo, S_dir_hi_lo, walkoff_angle_hi_lo, o_or_e_hi_lo, group_index_hi_lo, β0_hi_lo, β1_hi_lo, β2_hi_lo, β3_hi_lo
end

function calc_k_dir_ε_tensor_n_xyz(
    θ,
    ϕ,
    cr::NonlinearCrystal,
    λ::Unitful.Length=default_λ(cr),
    T::Unitful.Temperature=default_T(cr),
)
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"

    nx = refractive_index(cr.n_x_principal, λ, T)
    ny = refractive_index(cr.n_y_principal, λ, T)
    nz = refractive_index(cr.n_z_principal, λ, T)

    k_dir = calc_k_dir(θ, ϕ)
    ε_tensor = Diagonal(@SVector [nx^2, ny^2, nz^2])

    return k_dir, ε_tensor, (nx, ny, nz)
end

function calc_n_hi_lo_D_dir_hi_lo(k::AbstractVector{<:Real}, ε_tensor::AbstractMatrix{<:Real})
    P = I - k * k'
    ε_proj = P * ε_tensor * P

    # Use special native Julia eigen function to allow differentiation with ForwardDiff.jl
    eigs = DifferentiableEigen.eigen(ε_proj)
    eigvals = eigs[1][1:2:end]
    eigvecs = reshape(eigs[2][1:2:end], size(ε_proj))
    idx = sortperm(eigvals, rev=true)

    n_hi_lo = (sqrt(abs(eigvals[idx[1]])), sqrt(abs(eigvals[idx[2]])))
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

function assign_e_o(E_dir_hi_lo::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}}; angle_tol=0.2u"°")
    first_e = isapprox(abs(acos(clamp(dot([0, 0, 1], E_dir_hi_lo[1]), -1, 1))), 90.0u"°", atol=ustrip(u"rad", angle_tol))
    second_e = isapprox(abs(acos(clamp(dot([0, 0, 1], E_dir_hi_lo[2]), -1, 1))), 90.0u"°", atol=ustrip(u"rad", angle_tol))

    o_or_e_hi_lo = [nothing, nothing]
    if first_e && second_e
        o_or_e_hi_lo = [:o, :o]
    elseif second_e
        o_or_e_hi_lo = [:e, :o]
    elseif first_e
        o_or_e_hi_lo = [:o, :e]
    end
    return o_or_e_hi_lo
end

# β0 = k
function calc_β0_hi_lo(θ, ϕ, cr::NonlinearCrystal, λ::Unitful.Length=default_λ(cr), T::Unitful.Temperature=default_T(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    n_hi_lo = refraction_data_hi_lo(θ, ϕ, cr, λ, T; n_hi_lo_only=true)
    return Tuple([(2π / λ * n) |> u"m^-1" for n in n_hi_lo])
end
calc_β0_hi_lo(θ, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β0_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β1 = ∂k/∂ω
function calc_β1_hi_lo(θ, ϕ, cr::NonlinearCrystal, λ::Unitful.Length=default_λ(cr), T::Unitful.Temperature=default_T(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    ω_in = 2π * c_0 / λ

    sn = ω_in * 1e-3
    return (calc_β0_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in + sn), T) .- calc_β0_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in - sn), T)) ./ (2 * sn)

    # fun = ω -> [ustrip(d) for d in calc_β0_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)]
    # return Tuple(
    #     ForwardDiff.derivative(
    #         fun,
    #         ustrip(ω_in |> u"s^-1")
    #     ) * 1u"s/m")
end
calc_β1_hi_lo(θ, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β1_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β2 = ∂²k/∂ω²
function calc_β2_hi_lo(θ, ϕ, cr::NonlinearCrystal, λ::Unitful.Length=default_λ(cr), T::Unitful.Temperature=default_T(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    ω_in = 2π * c_0 / λ

    sn = ω_in * 1e-3
    (calc_β1_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in + sn), T) .- calc_β1_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in - sn), T)) ./ (2 * sn)

    # fun = ω -> [ustrip(d) for d in calc_β1_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)]
    # return Tuple(
    #     ForwardDiff.derivative(
    #         fun,
    #         ustrip(ω_in |> u"s^-1")
    #     ) * 1u"s^2/m")
end
calc_β2_hi_lo(θ, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β2_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

# β3 = ∂³k/∂ω³
function calc_β3_hi_lo(θ, ϕ, cr::NonlinearCrystal, λ::Unitful.Length=default_λ(cr), T::Unitful.Temperature=default_T(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    ω_in = 2π * c_0 / λ

    sn = ω_in * 1e-3
    (calc_β2_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in + sn), T) .- calc_β2_hi_lo(θ, ϕ, cr, 2π * c_0 / (ω_in - sn), T)) ./ (2 * sn)

    # fun = ω -> [ustrip(d) for d in calc_β2_hi_lo(θ, ϕ, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)]
    # return Tuple(
    #     ForwardDiff.derivative(
    #         fun,
    #         ustrip(ω_in |> u"s^-1")
    #     ) * 1u"s^3/m")
end
calc_β3_hi_lo(θ, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β3_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

function calc_group_index_hi_lo(θ, ϕ, cr::NonlinearCrystal, λ::Unitful.Length=default_λ(cr), T::Unitful.Temperature=default_T(cr))
    θ = θ |> u"rad"
    ϕ = ϕ |> u"rad"
    return Tuple([β1 * c_0 for β1 in calc_β1_hi_lo(θ, ϕ, cr, λ, T)])
end
calc_group_index_hi_lo(θ, cr::UnidirectionalCrystal, args...; kwargs...) = calc_group_index_hi_lo(θ, 0.0u"rad", cr, args...; kwargs...)

function plot_birefringent_refraction(
    θ,
    ϕ,
    cr::BidirectionalCrystal,
    λ::Unitful.Length=default_λ(cr),
    T::Unitful.Temperature=default_T(cr);
)
    k_dir, ε_tensor, (nx, ny, nz) = calc_k_dir_ε_tensor_n_xyz(θ, ϕ, cr, λ, T)
    n_hi_lo, D_dir_hi_lo, E_dir_hi_lo, S_dir_hi_lo, walkoff_angle_hi_lo = refraction_data_hi_lo(θ, ϕ, cr, λ, T)

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
    T::Union{AbstractVector{<:Unitful.Temperature},Unitful.Temperature}=[default_T(cr)],
    lambda_min=nothing,
    lambda_max=nothing,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")

    plot_refractiveindex!(cr.n_x_principal; n_sample_pts, T, lambda_min, lambda_max, label="n_x", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_y_principal; n_sample_pts, T, lambda_min, lambda_max, label="n_y", colormap=Reverse(:Greens))
    plot_refractiveindex!(cr.n_z_principal; n_sample_pts, T, lambda_min, lambda_max, label="n_z", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end

function plot_refractiveindex(
    cr::UnidirectionalCrystal;
    n_sample_pts=500,
    T::Union{AbstractVector{<:Unitful.Temperature},Unitful.Temperature}=[default_T(cr)],
    lambda_min=nothing,
    lambda_max=nothing,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")

    plot_refractiveindex!(cr.n_xy_principal; n_sample_pts, T, lambda_min, lambda_max, label="n_x/n_y (e)", colormap=Reverse(:Reds))
    plot_refractiveindex!(cr.n_z_principal; n_sample_pts, T, lambda_min, lambda_max, label="n_z (o)", colormap=Reverse(:Blues))
    Legend(f[1, 2], ax)
    DataInspector(ax)
    return f
end