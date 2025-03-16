export PhaseMatch, CollinearPhaseMatch, phase_match_wavelengths, find_phasematches, find_nearest_phase_match, delta_k, plot_delta_k_map

abstract type PhaseMatch end

struct CollinearPhaseMatch{LT,TT,OT,AT,IT,WT,DT,ET,ST,RT,GT,B2T,B3T,FT,CT}
    lambda_r1_r2_b::Vector{<:LT}
    T::TT
    hi_or_lo_r1_r2_b::Vector{Symbol}
    o_or_e_r1_r2_b::OT
    theta_pm::AT
    phi_pm::IT
    walkoff_angle_r1_r2_b::Vector{<:WT}
    D_dir_r1_r2_b::Vector{<:DT}
    E_dir_r1_r2_b::Vector{<:ET}
    S_dir_r1_r2_b::Vector{<:ST}
    refractive_index_r1_r2_b::Vector{<:RT}
    group_index_r1_r2_b::Vector{<:GT}
    beta_2_r1_r2_b::Vector{<:B2T}
    beta_3_r1_r2_b::Vector{<:B3T}
    d_eff::FT
    cr::CT
end

function CollinearPhaseMatch(
    cr::NonlinearCrystal,
    lambda_r1_r2_b::AbstractVector{<:Unitful.Length},
    T::Unitful.Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{<:Symbol},
    theta_pm,
    phi_pm,
)
    theta_pm = theta_pm |> u"rad"
    phi_pm = phi_pm |> u"rad"

    data_hi_lo_r1 = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[1], T)
    data_hi_lo_r2 = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[2], T)
    data_hi_lo_b = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[3], T)
    data_hi_lo = [data_hi_lo_r1, data_hi_lo_r2, data_hi_lo_b]

    select_hi_lo = data -> [(hi_or_lo_r1_r2_b[i] == :hi ? data[i][1] : data[i][2]) for i in eachindex(lambda_r1_r2_b)]
    refractive_index_r1_r2_b = select_hi_lo([d[1] for d in data_hi_lo])
    D_dir_r1_r2_b = select_hi_lo([d[2] for d in data_hi_lo])
    E_dir_r1_r2_b = select_hi_lo([d[3] for d in data_hi_lo])
    S_dir_r1_r2_b = select_hi_lo([d[4] for d in data_hi_lo])
    walkoff_angle_r1_r2_b = select_hi_lo([d[5] for d in data_hi_lo])
    o_or_e_r1_r2_b = select_hi_lo([d[6] for d in data_hi_lo])
    group_index_r1_r2_b = select_hi_lo([d[7] for d in data_hi_lo])
    # beta_0_r1_r2_b = select_hi_lo([d[8] for d in data_hi_lo])
    # beta_1_r1_r2_b = select_hi_lo([d[9] for d in data_hi_lo])
    beta_2_r1_r2_b = select_hi_lo([d[10] for d in data_hi_lo])
    beta_3_r1_r2_b = select_hi_lo([d[11] for d in data_hi_lo])

    d_eff = calc_d_eff(cr, E_dir_r1_r2_b...)

    return CollinearPhaseMatch(
        lambda_r1_r2_b,
        T,
        hi_or_lo_r1_r2_b,
        o_or_e_r1_r2_b,
        theta_pm,
        phi_pm,
        walkoff_angle_r1_r2_b,
        D_dir_r1_r2_b,
        E_dir_r1_r2_b,
        S_dir_r1_r2_b,
        refractive_index_r1_r2_b,
        group_index_r1_r2_b,
        beta_2_r1_r2_b,
        beta_3_r1_r2_b,
        d_eff,
        cr,
    )
end

function Base.show(io::IO, cpm::CollinearPhaseMatch)
    digits = 3
    println(io, "Crystal: $(cpm.cr.metadata[:description])")
    if !all(isnothing.(cpm.o_or_e_r1_r2_b))
        println(io, "Wavelengths: λ_r1: $(round(u"nm", cpm.lambda_r1_r2_b[1]; digits)) ($(cpm.hi_or_lo_r1_r2_b[1])/$(cpm.o_or_e_r1_r2_b[1])) + λ_r2: $(round(u"nm", cpm.lambda_r1_r2_b[2]; digits)) ($(cpm.hi_or_lo_r1_r2_b[2])/$(cpm.o_or_e_r1_r2_b[2])) = λ_b: $(round(u"nm", cpm.lambda_r1_r2_b[3]; digits)) ($(cpm.hi_or_lo_r1_r2_b[3])/$(cpm.o_or_e_r1_r2_b[3]))")
    else
        println(io, "Wavelengths: λ_r1: $(round(u"nm", cpm.lambda_r1_r2_b[1]; digits)) ($(cpm.hi_or_lo_r1_r2_b[1])) + λ_r2: $(round(u"nm", cpm.lambda_r1_r2_b[2]; digits)) ($(cpm.hi_or_lo_r1_r2_b[2])) = λ_b: $(round(u"nm", cpm.lambda_r1_r2_b[3]; digits)) ($(cpm.hi_or_lo_r1_r2_b[3]))")
    end
    println(io, "Temperature: $(cpm.T |> u"K") ($(float(cpm.T |> u"°C")))")
    println(io, "θ: $(round(u"°", cpm.theta_pm |> u"°"; digits)), ϕ: $(round(u"°", cpm.phi_pm |> u"°"; digits))")
    println(io, "Walkoff: $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[1] |> u"mrad"; digits)), $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[2] |> u"mrad"; digits)), $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[3] |> u"mrad"; digits))")
    println(io, "Refractive index: $(round(cpm.refractive_index_r1_r2_b[1]; digits)), $(round(cpm.refractive_index_r1_r2_b[2]; digits)), $(round(cpm.refractive_index_r1_r2_b[3]; digits))")
    println(io, "Group index: $(round(cpm.group_index_r1_r2_b[1]; digits)), $(round(cpm.group_index_r1_r2_b[2]; digits)), $(round(cpm.group_index_r1_r2_b[3]; digits))")
    println(io, "GDD: $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[1]; digits)), $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[2]; digits)), $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[3]; digits))")
    println(io, "TOD: $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[1]; digits)), $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[2]; digits)), $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[3]; digits))")
    println(io, "d_eff: $(round(u"pm/V", cpm.d_eff; digits))")
    println(io, "E_dir: $(round.(cpm.E_dir_r1_r2_b[1]; digits)), $(round.(cpm.E_dir_r1_r2_b[2]; digits)), $(round.(cpm.E_dir_r1_r2_b[3]; digits))")
end

function phase_match_wavelengths(;
    lambda_r1::Union{Nothing,Unitful.Length}=nothing,
    lambda_r2::Union{Nothing,Unitful.Length}=nothing,
    lambda_b::Union{Nothing,Unitful.Length}=nothing,
)
    if isnothing(lambda_r1)
        @assert (!isnothing(lambda_r2) && !isnothing(lambda_b))
        lambda_r1 = 1 / (1 / lambda_b - 1 / lambda_r2)
    end

    if isnothing(lambda_r2)
        @assert (!isnothing(lambda_r1) && !isnothing(lambda_b))
        lambda_r2 = 1 / (1 / lambda_b - 1 / lambda_r1)
    end

    if isnothing(lambda_b)
        @assert (!isnothing(lambda_r1) && !isnothing(lambda_r2))
        lambda_b = 1 / (1 / lambda_r1 + 1 / lambda_r2)
    end

    @assert 1 / lambda_r1 + 1 / lambda_r2 ≈ 1 / lambda_b

    return lambda_r1, lambda_r2, lambda_b
end

function delta_k(
    θ_pm,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::UnidirectionalCrystal;
    kwargs...
)
    return delta_k(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; kwargs...)
end

function delta_k(
    θ_pm,
    ϕ_pm,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Unitful.Length}=nothing,
    lambda_r2::Union{Nothing,Unitful.Length}=nothing,
    lambda_b::Union{Nothing,Unitful.Length}=nothing,
    T::Unitful.Temperature=cr.n_x_principal.T_ref
)
    θ_pm = θ_pm |> u"rad"
    ϕ_pm = ϕ_pm |> u"rad"
    lambda_r1, lambda_r2, lambda_b = phase_match_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_r1_r2_b])

    n_r1 = hi_or_lo_r1_r2_b[1] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, T; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, T; n_hi_lo_only=true)[2]
    n_r2 = hi_or_lo_r1_r2_b[2] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, T; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, T; n_hi_lo_only=true)[2]
    n_b = hi_or_lo_r1_r2_b[3] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, T; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, T; n_hi_lo_only=true)[2]

    return 2π * (n_r1 / lambda_r1 + n_r2 / lambda_r2 - n_b / lambda_b)
end

function delta_k(cpm::CollinearPhaseMatch)
    return delta_k(cpm.theta_pm, cpm.phi_pm, cpm.hi_or_lo_r1_r2_b, cpm.cr; lambda_r1=cpm.lambda_r1_r2_b[1], lambda_r2=cpm.lambda_r1_r2_b[2], lambda_b=cpm.lambda_r1_r2_b[3], cpm.T)
end

function find_phasematches(
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Unitful.Length}=nothing,
    lambda_r2::Union{Nothing,Unitful.Length}=nothing,
    lambda_b::Union{Nothing,Unitful.Length}=nothing,
    T::Unitful.Temperature=default_T(cr),
    theta_fixed=nothing,
    phi_fixed=nothing,
    ngrid=500,
    tol=1e-6
)
    lambda_r1, lambda_r2, lambda_b = phase_match_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    lambda_r1_r2_b = [lambda_r1, lambda_r2, lambda_b]

    phasematches = CollinearPhaseMatch[]

    if theta_fixed !== nothing && isnothing(phi_fixed)
        # Global search over ϕ
        all_ϕ = range(0, 2π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T) for ϕ in all_ϕ]

        # Find approximate zero-crossings
        for i in 1:(length(all_ϕ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0  # zero crossing detected
                # Refine solution using local optimization
                ϕ_sol = find_zero(ϕ -> delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T),
                    (all_ϕ[i], all_ϕ[i+1]), Bisection(), atol=tol)
                push!(phasematches, CollinearPhaseMatch(cr, lambda_r1_r2_b, T, hi_or_lo_r1_r2_b, theta_fixed, ϕ_sol))
            end
        end

    elseif isnothing(theta_fixed) && !isnothing(phi_fixed)
        all_θ = range(0, π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T) for θ in all_θ]

        # Find approximate zero-crossings
        for i in 1:(length(all_θ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0 # zero crossing detected
                # Refine solution using local optimization
                θ_sol = find_zero(θ -> delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T),
                    (all_θ[i], all_θ[i+1]), Bisection(), atol=tol)
                push!(phasematches, CollinearPhaseMatch(cr, lambda_r1_r2_b, T, hi_or_lo_r1_r2_b, θ_sol, phi_fixed))
            end
        end
    else
        error("You must provide either theta_fixed or phi_fixed.")
    end
    return phasematches
end

function find_nearest_phase_match(
    θ_target,
    ϕ_target,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Unitful.Length}=nothing,
    lambda_r2::Union{Nothing,Unitful.Length}=nothing,
    lambda_b::Union{Nothing,Unitful.Length}=nothing,
    T::Unitful.Temperature=default_T(cr),
    ngrid=500,
    tol=1e-6
)
    θ_target = θ_target |> u"rad"
    ϕ_target = ϕ_target |> u"rad"

    all_pm_θ_fixed = find_phasematches(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T, ngrid, tol, theta_fixed=θ_target)
    all_pm_ϕ_fixed = find_phasematches(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T, ngrid, tol, phi_fixed=ϕ_target)

    all_pm_candidates = [[pm for pm in all_pm_θ_fixed]; [pm for pm in all_pm_ϕ_fixed]]
    isempty(all_pm_candidates) && return nothing

    all_θ_pm_candidates = [pm.theta_pm for pm in all_pm_candidates]
    all_ϕ_pm_candidates = [pm.phi_pm for pm in all_pm_candidates]

    target_vec = angles_to_vector(θ_target, ϕ_target)
    pm_candidate_vecs = angles_to_vector.(all_θ_pm_candidates, all_ϕ_pm_candidates)
    i_nearest = findmax([dot(pm_vec, target_vec) for pm_vec in pm_candidate_vecs])[2]

    return all_pm_candidates[i_nearest]
end

function plot_delta_k_map(
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Unitful.Length}=nothing,
    lambda_r2::Union{Nothing,Unitful.Length}=nothing,
    lambda_b::Union{Nothing,Unitful.Length}=nothing,
    T::Unitful.Temperature=default_T(cr),
    n_points::Integer=100,
    axis_length::Real=1.5,
    digits::Integer=3,
    plot_type::Symbol=:polar,
    theta_fixed=nothing,
    phi_fixed=nothing,
    show_coordinates::Bool=true,
    show_optical_axes::Bool=true,
)
    # Color settings
    col_coordinates = (color=5, colormap=:Greys, colorrange=(1, 10))
    col_contour = (color=1, colormap=:Set1_9, colorrange=(1, 9))
    col_r1 = (color=9, colormap=:vik10, colorrange=(1, 10))
    col_r2 = (color=8, colormap=:vik10, colorrange=(1, 10))
    col_b = (color=3, colormap=:vik10, colorrange=(1, 10))
    col_heatmap = (colormap=:vik,)

    # Prepare data
    lambda_r1_r2_b = phase_match_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    if isnothing(T)
        T = isa(cr, UnidirectionalCrystal) ? cr.n_o_principal.T_ref : cr.n_x_principal.T_ref
    end

    θ_range, ϕ_range, all_delta_k = compute_delta_k_grid(cr, hi_or_lo_r1_r2_b, lambda_r1_r2_b..., T, n_points)
    scale_limit = maximum(abs.(all_delta_k))

    if !isnothing(phi_fixed) || !isnothing(theta_fixed)
        phasematches = find_phasematches(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T, theta_fixed, phi_fixed)
    else
        phasematches = []
    end

    # Generate contours
    cpl = Makie.Contours.contours(ustrip.(θ_range), ustrip.(ϕ_range), all_delta_k, [0.0])

    # Compute optical axes
    oa = [optical_axis_angle(cr, λ, T) for λ in lambda_r1_r2_b]

    # Call appropriate plotting subfunction
    if plot_type == :polar
        return plot_polar_mode(
            θ_range,
            ϕ_range,
            all_delta_k,
            scale_limit,
            cr,
            T,
            cpl,
            phasematches,
            oa,
            lambda_r1_r2_b,
            hi_or_lo_r1_r2_b,
            digits,
            show_coordinates,
            show_optical_axes,
            col_coordinates,
            col_contour,
            col_r1,
            col_r2,
            col_b,
            col_heatmap,
        )
    elseif plot_type == :sphere
        return plot_sphere_mode(
            θ_range,
            ϕ_range,
            all_delta_k,
            scale_limit,
            cr,
            T,
            cpl,
            phasematches,
            oa,
            lambda_r1_r2_b,
            hi_or_lo_r1_r2_b,
            digits,
            axis_length,
            show_coordinates,
            show_optical_axes,
            col_coordinates,
            col_contour,
            col_r1,
            col_r2,
            col_b,
            col_heatmap,
        )
    else
        error("Plot type $(plot_type) unknown.")
    end
end

function compute_delta_k_grid(
    cr::NonlinearCrystal,
    hi_or_lo_r1_r2_b,
    lambda_r1,
    lambda_r2,
    lambda_b,
    T,
    n_points,
)
    θ_range = LinRange(0, π, n_points) * u"rad"
    ϕ_range = LinRange(0, 2π, 2 * n_points) * u"rad"

    # TODO: Use symmetry information from cr.metadata[:pointgroup] for speedup
    if isa(cr, UnidirectionalCrystal)
        # No ϕ dependence, compute only one value and repeat
        all_delta_k = [
            delta_k(θ, 0.0u"°", hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T)
            for θ in θ_range
        ]
        all_delta_k = repeat(reshape(all_delta_k, :, 1), 1, length(ϕ_range))
    else
        all_delta_k = [delta_k(θ, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, T) for θ in θ_range, ϕ in ϕ_range]
    end

    return θ_range, ϕ_range, ustrip.(u"µm^-1", all_delta_k)
end


function plot_polar_mode(
    θ_range::AbstractVector,
    ϕ_range::AbstractVector,
    all_delta_k::AbstractMatrix{<:Real},
    scale_limit::Real,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    cpl,
    phasematches::AbstractVector,
    oa::AbstractVector,
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length},
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    digits::Integer,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    col_coordinates,
    col_contour,
    col_r1,
    col_r2,
    col_b,
    col_heatmap,
)
    f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel="ϕ",
        ylabel="θ",
        title="$(round(u"nm", lambda_r1_r2_b[1]; digits)) (λ_r1, $(hi_or_lo_r1_r2_b[1])) + $(round(u"nm", lambda_r1_r2_b[2]; digits)) (λ_r2, $(hi_or_lo_r1_r2_b[2])) = $(round(u"nm", lambda_r1_r2_b[3]; digits)) (λ_b, $(hi_or_lo_r1_r2_b[3]))"
    )

    heatmap!(ax, ϕ_range .|> u"°", θ_range .|> u"°", all_delta_k'; col_heatmap..., colorrange=(-scale_limit, scale_limit), inspectable=false)

    plot_polar_contours!(ax, cpl, col_contour, cr, T, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_polar_phasematches!(ax, phasematches)

    if show_coordinates
        plot_polar_coordinate_markers!(ax, col_coordinates)
    end

    if show_optical_axes
        plot_polar_optical_axes!(ax, oa; col_r1, col_r2, col_b)
    end

    # Legend(f[1, 2], ax)
    Colorbar(f[1, 2]; col_heatmap..., limits=(-scale_limit, scale_limit), label="Δk in µm⁻¹")
    DataInspector(ax)

    return f
end

function plot_phasematch_label(pm::CollinearPhaseMatch)
    buf = IOBuffer()
    print(buf, pm)
    return String(take!(buf))
end

function polar_plot_phasematch_label(
    plot,
    index,
    position,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length}
)
    θ_pm = position[2] * u"°"
    ϕ_pm = position[1] * u"°"
    pm = find_nearest_phase_match(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; T, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
    return plot_phasematch_label(pm)
end

function plot_polar_contours!(
    ax::Axis,
    cpl,
    col_contour,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length};
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θi_cl = [v[1] for v in cl.vertices] .|> u"°"
        ϕi_cl = [v[2] for v in cl.vertices] .|> u"°"

        label_fun = (plot, index, position) -> polar_plot_phasematch_label(plot, index, position, cr, T, hi_or_lo_r1_r2_b, lambda_r1_r2_b)

        lines!(
            ax, ϕi_cl, θi_cl;
            color=:black, linewidth=2, fxaa=true, inspectable=true,
            inspector_label=label_fun,
            label=(i == 1 ? "Phase matches" : nothing), col_contour...
        )
    end
end

function plot_polar_phasematches!(ax::Axis, phasematches::AbstractVector)
    for (i, pm) in enumerate(phasematches)
        scatter!(
            ax, pm.phi_pm |> u"°", pm.theta_pm |> u"°";
            markersize=10,
            label="Probed phase match $(i)", inspectable=true,
            inspector_label=(plot, index, position) -> ("Probed phase match $(i)\n" * plot_phasematch_label(pm))
        )
    end
end

function plot_polar_coordinate_markers!(ax::Axis, col_coordinates)
    digits = 3
    all_coordinate_labels = [
        (0.0, 0.0, "Z", (:left, :bottom)),
        (0.0, π / 2, "Z", (:left, :bottom)),
        (0.0, π, "Z", (:left, :bottom)),
        (0.0, 3π / 2, "Z", (:left, :bottom)),
        (0.0, 2π, "Z", (:right, :bottom)),
        (π / 2, 0.0, "X", (:left, :bottom)),
        (π / 2, π / 2, "Y", (:left, :bottom)),
        (π / 2, π, "X", (:left, :bottom)),
        (π / 2, 3π / 2, "Y", (:left, :bottom)),
        (π / 2, 2π, "X", (:right, :bottom)),
        (π, 0.0, "Z", (:left, :top)),
        (π, π / 2, "Z", (:right, :top)),
        (π, π, "Z", (:right, :top)),
        (π, 3π / 2, "Z", (:right, :top)),
        (π, 2π, "Z", (:right, :top)),
    ]
    for (θ, ϕ, lab, a) in all_coordinate_labels
        θ = θ |> u"°"
        ϕ = ϕ |> u"°"
        coord_label = "$(lab) axis\nθ = $(round(u"°", θ; digits)) ($(round(u"rad", θ; digits)))\nϕ = $(round(u"°", ϕ; digits)) ($(round(u"rad", ϕ; digits)))"
        scatter!(ax, ϕ, θ; col_coordinates..., inspectable=true, inspector_label=(plot, index, position) -> coord_label)
        text!(ax, ϕ, θ; text=lab, col_coordinates..., align=a)
    end
end

function plot_polar_optical_axes!(ax::Axis, oa; col_r1, col_r2, col_b)
    digits = 3
    colors = [col_r1, col_r2, col_b]
    labels = ["λ_r1", "λ_r2", "λ_b"]

    for (angle, col, lab) in zip(oa, colors, labels)
        axis_label = "Optical axis $(lab)\n|θ| = $(round(u"°", angle; digits)) ($(round(u"rad", angle; digits)))\nϕ = $(0.0u"°") ($(0.0u"rad"))"
        scatter!(
            ax,
            [0, 0, π, π, 2π, 2π] .|> u"°",
            [angle, π - angle, angle, π - angle, angle, π - angle] .|> u"°";
            col..., markersize=8, label=lab,
            inspectable=true, inspector_label=(plot, index, position) -> axis_label,
        )
    end
end


function plot_sphere_mode(
    θ_range,
    ϕ_range,
    all_delta_k::AbstractMatrix{<:Real},
    scale_limit::Real,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    cpl,
    phasematches::AbstractVector,
    oa::AbstractVector,
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length},
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    digits::Integer,
    axis_length::Real,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    col_coordinates,
    col_contour,
    col_r1,
    col_r2,
    col_b,
    col_heatmap,
)
    f = Figure()
    ax = Axis3(
        f[1, 1],
        azimuth=0.1π, elevation=0.05π,
        aspect=:data, viewmode=:fit,
        title="$(round(u"nm", lambda_r1_r2_b[1]; digits)) (λ_r1, $(hi_or_lo_r1_r2_b[1])) + $(round(u"nm", lambda_r1_r2_b[2]; digits)) (λ_r2, $(hi_or_lo_r1_r2_b[2])) = $(round(u"nm", lambda_r1_r2_b[3]; digits)) (λ_b, $(hi_or_lo_r1_r2_b[3]))"
    )
    hidedecorations!(ax, grid=false)

    plot_sphere_mesh!(ax, θ_range, ϕ_range, all_delta_k, scale_limit, col_heatmap)
    plot_sphere_contours!(ax, cpl, col_contour, cr, T, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_sphere_phasematches!(ax, phasematches, axis_length; col_r1, col_r2, col_b)

    if show_coordinates
        plot_sphere_coordinates!(ax, axis_length, col_coordinates)
    end

    if show_optical_axes
        plot_sphere_optical_axes!(ax, oa, axis_length; col_r1, col_r2, col_b)
    end

    # Legend(f[1, 2], ax)
    Colorbar(
        f[1, 2]; col_heatmap...,
        limits=(-scale_limit, scale_limit),
        label="Δk in µm⁻¹"
    )
    DataInspector(ax)

    return f
end

function plot_sphere_mesh!(ax::Axis3, θ_range, ϕ_range, all_delta_k, scale_limit, col_heatmap)
    points = [Point3f(angles_to_vector(θ, ϕ)...) for θ in θ_range, ϕ in ϕ_range]
    faces = decompose(QuadFace{GLIndex}, Tessellation(Rect(0, 0, 1, 1), size(all_delta_k)))
    normals = normalize.(vec(points))
    gb_mesh = GeometryBasics.Mesh(vec(points), faces; normal=normals)
    mesh!(ax, gb_mesh, color=all_delta_k[:]; col_heatmap..., colorrange=(-scale_limit, scale_limit), inspectable=false)
end

function sphere_plot_phasematch_label(
    plot,
    index,
    position,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length},
)
    θ_pm, ϕ_pm = vector_to_angles(position)
    pm = find_nearest_phase_match(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; T, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
    return plot_phasematch_label(pm)
end

function plot_sphere_contours!(
    ax::Axis3,
    cpl,
    col_contour,
    cr::NonlinearCrystal,
    T::Unitful.Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Unitful.Length,Unitful.Length,Unitful.Length}
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θ = [v[1] for v in cl.vertices]
        ϕ = [v[2] for v in cl.vertices]

        label_fun = (plot, index, position) -> sphere_plot_phasematch_label(plot, index, position, cr, T, hi_or_lo_r1_r2_b, lambda_r1_r2_b)

        points = [Point3f(angles_to_vector(θi, ϕi)...) for (θi, ϕi) in zip(θ, ϕ)]
        lines!(ax, points; col_contour..., linewidth=4, label=(i == 1 ? "Phase matches" : nothing), inspectable=true, inspector_label=label_fun,
        )
    end
end

function plot_sphere_phasematches!(
    ax::Axis3,
    phasematches,
    axis_length;
    col_r1,
    col_r2,
    col_b
)
    fac_r1 = [1 + 2 * (axis_length - 1) / 3, axis_length]
    fac_r2 = [1 + (axis_length - 1) / 3, 1 + 2 * (axis_length - 1) / 3]
    fac_b = [1, 1 + (axis_length - 1) / 3]

    len_E_dir = (axis_length - 1) / 3

    lambda_labels = ["λ_r1", "λ_r2", "λ_b"]
    factors = [fac_r1, fac_r2, fac_b]
    colors = [col_r1, col_r2, col_b]

    for (i, pm) in enumerate(phasematches)
        θi_pm = pm.theta_pm
        ϕi_pm = pm.phi_pm
        E_dirs = pm.E_dir_r1_r2_b

        ri_pm = Point3f(angles_to_vector(θi_pm, ϕi_pm)...)

        scatter!(ax, axis_length .* ri_pm; markersize=15, label="Probed phase match $(i)", inspectable=true,
            inspector_label=(plot, index, position) -> ("Probed phase match $(i)\n" * plot_phasematch_label(pm)))

        for (i, (lab, fac, E_dir, col)) in enumerate(zip(lambda_labels, factors, E_dirs, colors))
            label_fun = (plot, index, position) -> ("$(lab)")
            lines!(
                ax,
                [fac[1] .* ri_pm, fac[2] .* ri_pm];
                col..., linewidth=4, inspectable=true, inspector_label=label_fun
            )

            scale = sum(fac) / 2
            lines!(
                ax,
                scale .* ri_pm[1] .+ len_E_dir .* E_dir[1] .* [-1, +1],
                scale .* ri_pm[2] .+ len_E_dir .* E_dir[2] .* [-1, +1],
                scale .* ri_pm[3] .+ len_E_dir .* E_dir[3] .* [-1, +1];
                col..., linewidth=4, inspectable=true, inspector_label=label_fun
            )
        end
    end
end



function plot_sphere_coordinates!(ax::Axis3, axis_length, col_coordinates)
    digits = 3
    coords = [
        ([-1, 1], [0, 0], [0, 0], "X"),
        ([0, 0], [-1, 1], [0, 0], "Y"),
        ([0, 0], [0, 0], [-1, 1], "Z")
    ]
    for (x, y, z, label) in coords
        θ, ϕ = vector_to_angles([x[2], y[2], z[2]])
        θ = θ |> u"°"
        ϕ = ϕ |> u"°"
        coord_label = "$(label) axis\nθ = $(round(u"°", θ; digits)) ($(round(u"rad", θ; digits)))\nϕ = $(round(u"°", ϕ; digits)) ($(round(u"rad", ϕ; digits)))"
        lines!(ax, axis_length .* x, axis_length .* y, axis_length .* z; col_coordinates..., inspectable=true, inspector_label=(plot, index, position) -> coord_label, linewidth=3)
        text!(ax, axis_length .* x[2], axis_length .* y[2], axis_length .* z[2]; text=label, col_coordinates...)
        text!(ax, axis_length .* x[1], axis_length .* y[1], axis_length .* z[1]; text=label, col_coordinates...)
    end
end

function plot_sphere_optical_axes!(ax::Axis3, oa, axis_length; col_r1, col_r2, col_b)
    digits = 3
    colors = [col_r1, col_r2, col_b]
    labels = ["λ_r1", "λ_r2", "λ_b"]

    for (angle, col, lab) in zip(oa, colors, labels)
        axis_label = "Optical axis $(lab)\n|θ| = $(round(u"°", angle; digits)) ($(round(u"rad", angle; digits)))\nϕ = $(0.0u"°") ($(0.0u"rad"))"
        dir_pos = axis_length * Point3f(sin(angle), 0, cos(angle))
        dir_neg = -dir_pos
        lines!(ax, [dir_neg, dir_pos]; col..., linewidth=4, label="Optical axes " * lab, inspectable=true, inspector_label=(plot, index, position) -> axis_label)

        dir_pos_neg = axis_length * Point3f(sin(-angle), 0, cos(-angle))
        dir_neg_neg = -dir_pos_neg
        lines!(ax, [dir_neg_neg, dir_pos_neg]; col..., linewidth=4, inspectable=true, inspector_label=(plot, index, position) -> axis_label)
    end
end