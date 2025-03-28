export calc_noncritical_pm_lines, plot_noncritical_pms, plot_delta_k_map

function noncritical_pm_label(
    plot,
    index,
    position,
    θ_pm,
    ϕ_pm,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r_switch::Symbol,
)
    lambda_b = position[1] * u"µm"
    lambda_r1 = lambda_r2 = nothing
    if lambda_r_switch == :lambda_r1
        lambda_r1 = position[2] * u"µm"
    elseif lambda_r_switch == :lambda_r2
        lambda_r2 = position[2] * u"µm"
    end
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b)
    return isnothing(pm) ? "" : plot_pm_label(pm)
end

function calc_raw_noncritical_pm_lines(
    principal_axis::Symbol,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
)
    @assert principal_axis in [:X, :Y, :Z]

    isnothing(lambda_b_min) && (lambda_b_min = valid_lambda_range(cr)[1])
    isnothing(lambda_b_max) && (lambda_b_max = valid_lambda_range(cr)[2] / 2)
    isnothing(lambda_r12_min) && (lambda_r12_min = valid_lambda_range(cr)[1])
    isnothing(lambda_r12_max) && (lambda_r12_max = valid_lambda_range(cr)[2])
    range_lambda_b = LinRange(lambda_b_min, lambda_b_max, ngrid) .|> u"µm"
    range_lambda_r12 = LinRange(lambda_r12_min, lambda_r12_max, ngrid) .|> u"µm"

    all_delta_k = zeros(length(range_lambda_b), length(range_lambda_r12))
    for i in CartesianIndices(all_delta_k)
        (i_b, i_r12) = (i[1], i[2])
        lambda_b = range_lambda_b[i_b]
        lambda_r12 = range_lambda_r12[i_r12]
        lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_b, lambda_r1=lambda_r12)
        if is_lambda_valid(lambda_r2, cr)
            all_delta_k[i_b, i_r12] = ustrip(u"m^-1", delta_k(axes_to_θ_ϕ(principal_axis)[1]..., hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b))
        else
            all_delta_k[i_b, i_r12] = NaN
        end
    end

    # Extract Δk=0 isolines using Makie contour function
    cpl = Makie.Contours.contour(range_lambda_b, range_lambda_r12, all_delta_k, 0.0)
    return cpl # Makie.jl contour object with list of `lines` containing list of `vertices`
end

function split_and_interpolate(cb::AbstractVector, cr::AbstractVector, sign_switches, fractions, segment_signs)
    segments_cb = []
    segments_cr = []

    start_idx = 1
    cb_intersections = []
    cr_intersections = []
    for (switch_idx, frac) in zip(sign_switches, fractions)
        # Extract current segment (from start_idx to switch_idx)
        cb_seg = cb[start_idx:switch_idx]
        cr_seg = cr[start_idx:switch_idx]

        # Interpolate at switch point
        cb_interp = cb[switch_idx] * (1 - frac) + cb[switch_idx+1] * frac
        cr_interp = cr[switch_idx] * (1 - frac) + cr[switch_idx+1] * frac

        # Add last point from previous segment at start of the current
        !isempty(cb_intersections) && prepend!(cb_seg, cb_intersections[end])
        !isempty(cb_intersections) && prepend!(cr_seg, cr_intersections[end])

        # Append interpolated point
        push!(cb_seg, cb_interp)
        push!(cr_seg, cr_interp)

        # Save segment
        push!(segments_cb, cb_seg)
        push!(segments_cr, cr_seg)

        # Start the next segment with the last value
        push!(cb_intersections, cb_interp)
        push!(cr_intersections, cr_interp)

        # Update start index
        start_idx = switch_idx + 1
    end

    # Last segment
    cb_seg = cb[start_idx:end]
    cr_seg = cr[start_idx:end]
    !isempty(cb_intersections) && prepend!(cb_seg, cb_intersections[end])
    !isempty(cr_intersections) && prepend!(cr_seg, cr_intersections[end])
    push!(segments_cb, cb_seg)
    push!(segments_cr, cr_seg)

    return segments_cb, segments_cr, cb_intersections, cr_intersections, segment_signs
end

function calc_noncritical_pm_lines(
    principal_axis::Symbol,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
)
    cpl = calc_raw_noncritical_pm_lines(
        principal_axis,
        hi_or_lo_r1_r2_b,
        cr;
        lambda_b_min,
        lambda_b_max,
        lambda_r12_min,
        lambda_r12_max,
        temp,
        ngrid
    )

    # If both red waves are polarized equally (type 1 phasematching), the plot lines show a helpful phase match symmetry around the line cont_r_raw = 2 * cont_b_raw
    is_type_one = (hi_or_lo_r1_r2_b[1] == hi_or_lo_r1_r2_b[2])

    all_segments_r = []
    all_segments_b = []
    all_cb_intersections = []
    all_cr_intersections = []
    # Iterate over all identified raw contour lines and postprocess/refine them
    for cl in cpl.lines
        cont_b_raw = [v[1] for v in cl.vertices]
        cont_r_raw = [v[2] for v in cl.vertices]

        # Kick out all NaNs and split contours into segments at all NaNs
        for (cont_b, cont_r) in zip(split_on_nan(cont_b_raw, cont_r_raw)...)
            shear_cr = cont_r .- 2 * cont_b # All SHG points are now on the y = 0 axis
            segments_cb, segments_cr, cb_intersections, cr_intersections, segment_signs = split_and_interpolate(cont_b, cont_r, sign_switch_fractions(shear_cr)...)

            if is_type_one
                # Throw away lower (r2) part, as it is symmetric anyway
                r1_segments = findall(x -> x > 0.0, segment_signs)
                segments_cb = segments_cb[r1_segments]
                segments_cr = segments_cr[r1_segments]
            else
                # Flip lower (r2) parts to upper (r1) half using phase match symmetry
                segments_cr = [segment_signs[i] < 0.0 ? 1 ./ (1 ./ segments_cb[i] - 1 ./ segments_cr[i]) : segments_cr[i] for i in eachindex(segments_cr)]
            end
            push!(all_segments_r, segments_cr)
            push!(all_segments_b, segments_cb)
            push!(all_cb_intersections, cb_intersections)
            push!(all_cr_intersections, cr_intersections)
        end
    end

    all_segments_b = reduce(vcat, all_segments_b; init=[])
    all_segments_r = reduce(vcat, all_segments_r; init=[])
    all_cb_intersections = reduce(vcat, all_cb_intersections; init=[])
    all_cr_intersections = reduce(vcat, all_cr_intersections; init=[])
    return all_segments_b, all_segments_r, all_cb_intersections, all_cr_intersections
end

function plot_single_noncritical_pm!(
    ax::Axis,
    principal_axis::Symbol,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=50,
    tol=1e-14u"nm^-1",
)
    is_type_1 = hi_or_lo_r1_r2_b[1] === hi_or_lo_r1_r2_b[2]

    all_segments_b, all_segments_r1, all_cb_intersections, all_cr_intersections = calc_noncritical_pm_lines(
        principal_axis,
        hi_or_lo_r1_r2_b,
        cr;
        lambda_b_min,
        lambda_b_max,
        lambda_r12_min,
        lambda_r12_max,
        temp,
        ngrid,
    )

    vline_lab = (plot, idx, pos) -> begin
        pos = pos * u"µm"
        marker_x[] = [ustrip(u"µm", pos[1])]
        marker_x_b[] = [pos[1]]
        is_r1 = pos[2] > 2 * pos[1]
        marker_y_r1[] = is_r1 ? [pos[2]] : [1 / (1 / pos[1] - 1 / pos[2])]
        marker_y_r2[] = is_r1 ? [1 / (1 / pos[1] - 1 / pos[2])] : [pos[2]]

        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_r1_r2_b, cr; lambda_r1=marker_y_r1[][1], lambda_b=pos[1], temp, principal_axis, ngrid, tol)
        return "$(is_type_1 ? "Type 1" : "Type 2") noncritical phasematch\n" * plot_pm_label(pm)
    end

    shg_lab = (plot, idx, pos) -> begin
        pos = pos * u"µm"
        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_r1_r2_b, cr; lambda_r1=pos[2], lambda_b=pos[1], temp, principal_axis, ngrid, tol)
        return "$(is_type_1 ? "Type 1" : "Type 2") noncritical phasematch SHG point\n" * plot_pm_label(pm)
    end

    vline_clear = (inspector, plot) -> begin
        marker_x[] = [NaN]
        marker_x_b[] = [NaN * u"µm"]
        marker_y_r1[] = [NaN * u"µm"]
        marker_y_r2[] = [NaN * u"µm"]
        return nothing
    end

    [lines!(ax, all_segments_b[i], all_segments_r1[i]; linestyle=hi_or_lo_r1_r2_b[1] === :hi ? :dash : :solid, color=COL_R1, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
    [lines!(ax, all_segments_b[i], 1 ./ (1 ./ all_segments_b[i] .- 1 ./ all_segments_r1[i]); linestyle=hi_or_lo_r1_r2_b[2] === :hi ? :dash : :solid, color=COL_R2, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
    !isempty(all_cb_intersections) && scatter!(ax, all_cb_intersections, all_cr_intersections; inspectable=true, inspector_label=shg_lab)

    marker_x = Observable([NaN]) # Workaround: vlines! is not yet compatible with unit observables
    marker_x_b = Observable([NaN * u"µm"])
    marker_y_r1 = Observable([NaN * u"µm"])
    marker_y_r2 = Observable([NaN * u"µm"])
    vlines!(ax, marker_x; inspectable=false, color=:gray, linewidth=2, alpha=0.5)
    scatter!(ax, marker_x_b, marker_y_r1; inspectable=false, color=COL_R1)
    scatter!(ax, marker_x_b, marker_y_r2; inspectable=false, color=COL_R2)
    return ax
end

function plot_noncritical_pms(
    principal_axis::Symbol,
    cr::NonlinearCrystal;
    hi_or_lo_r1_r2_b::Union{AbstractVector{Symbol},AbstractVector{<:AbstractVector{Symbol}},Nothing}=nothing,
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=50,
    tol=1e-14u"nm^-1",
    size::NTuple{2,Int}=(800, 600),
)

    f = Figure(; size)
    uc = Makie.UnitfulConversion(u"µm"; units_in_label=true)
    ax = Axis(
        f[1, 1],
        xlabel="λ_b",
        ylabel="λ_r12",
        title="$(cr.metadata[:description])\nNoncritical phasematches along positive $(principal_axis) axis, Temperature: $(float(temp |> u"K")) ($(float(temp |> u"°C")))", # , polarization directions: $(hi_or_lo_r1_r2_b) # TODO: Add lambda names and o/e
        dim1_conversion=uc,
        dim2_conversion=uc,
    )

    if isnothing(hi_or_lo_r1_r2_b)
        hi_or_lo_r1_r2_b = bool_permutations(:hi, :lo, 3)
        setdiff!(hi_or_lo_r1_r2_b, [[:hi, :lo, :lo], [:hi, :lo, :hi], [:lo, :lo, :lo], [:hi, :hi, :hi]]) # Prevent double plotting of Type 2 noncritical phasematches and unphysical phasematches 
    elseif typeof(hi_or_lo_r1_r2_b) <: AbstractVector{Symbol}
        hi_or_lo_r1_r2_b = [hi_or_lo_r1_r2_b]
    end

    for hl in hi_or_lo_r1_r2_b
        plot_single_noncritical_pm!(
            ax,
            principal_axis,
            hl,
            cr;
            lambda_b_min,
            lambda_b_max,
            lambda_r12_min,
            lambda_r12_max,
            temp,
            ngrid,
            tol,
        )
    end

    DataInspector(f)
    return f
end

function plot_delta_k_map(
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    pms::AbstractVector{<:PhaseMatch}=PhaseMatch[],
    n_points::Integer=100,
    axis_length::Real=1.5,
    digits::Integer=3,
    plot_type::Symbol=:polar,
    show_coordinates::Bool=true,
    show_optical_axes::Bool=true,
    size::NTuple{2,Int}=(800, 600),
)
    # Prepare data
    lambda_r1_r2_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    if isnothing(temp)
        temp = isa(cr, UnidirectionalCrystal) ? cr.n_o_principal.temp_ref : cr.n_x_principal.temp_ref
    end

    θ_range, ϕ_range, all_delta_k = compute_delta_k_grid(cr, hi_or_lo_r1_r2_b, lambda_r1_r2_b..., temp, n_points)
    scale_limit = maximum(abs.(all_delta_k))

    # Generate contours
    cpl = Makie.Contours.contours(ustrip.(θ_range), ustrip.(ϕ_range), all_delta_k, [0.0])

    # Compute optical axes
    oa = [optical_axis_angle(cr, λ, temp) for λ in lambda_r1_r2_b]

    # Call appropriate plotting subfunction
    if plot_type == :polar
        return plot_polar_mode(
            θ_range,
            ϕ_range,
            all_delta_k,
            scale_limit,
            cr,
            temp,
            cpl,
            pms,
            oa,
            lambda_r1_r2_b,
            hi_or_lo_r1_r2_b,
            digits,
            show_coordinates,
            show_optical_axes,
            size,
        )
    elseif plot_type == :sphere
        return plot_sphere_mode(
            θ_range,
            ϕ_range,
            all_delta_k,
            scale_limit,
            cr,
            temp,
            cpl,
            pms,
            oa,
            lambda_r1_r2_b,
            hi_or_lo_r1_r2_b,
            digits,
            axis_length,
            show_coordinates,
            show_optical_axes,
            size,
        )
    else
        error("Plot type $(plot_type) unknown.")
    end
end

plot_delta_k_map(args...; pms::PhaseMatch, kwargs...) = plot_delta_k_map(args...; pms=[pms], kwargs...)

function plot_delta_k_map(
    pm::PhaseMatch;
    n_points::Integer=100,
    axis_length::Real=1.5,
    digits::Integer=3,
    plot_type::Symbol=:polar,
    show_coordinates::Bool=true,
    show_optical_axes::Bool=true,
)
    return plot_delta_k_map(pm.hi_or_lo_r1_r2_b, pm.cr; lambda_r1=pm.lambda_r1_r2_b[1], lambda_r2=pm.lambda_r1_r2_b[2], lambda_b=pm.lambda_r1_r2_b[3], temp=pm.temp, pms=[pm], n_points, axis_length, digits, plot_type, show_coordinates, show_optical_axes)
end

function compute_delta_k_grid(
    cr::NonlinearCrystal,
    hi_or_lo_r1_r2_b,
    lambda_r1,
    lambda_r2,
    lambda_b,
    temp,
    n_points,
)
    θ_range = LinRange(0, π, n_points) * u"rad"
    ϕ_range = LinRange(0, 2π, 2 * n_points) * u"rad"

    # TODO: Use symmetry information from cr.metadata[:pointgroup] for speedup
    if isa(cr, UnidirectionalCrystal)
        # No ϕ dependence, compute only one value and repeat
        all_delta_k = [
            delta_k(θ, 0.0u"°", hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp)
            for θ in θ_range
        ]
        all_delta_k = repeat(reshape(all_delta_k, :, 1), 1, length(ϕ_range))
    else
        all_delta_k = [delta_k(θ, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in θ_range, ϕ in ϕ_range]
    end

    return θ_range, ϕ_range, ustrip.(u"µm^-1", all_delta_k)
end


function plot_polar_mode(
    θ_range::AbstractVector,
    ϕ_range::AbstractVector,
    all_delta_k::AbstractMatrix{<:Real},
    scale_limit::Real,
    cr::NonlinearCrystal,
    temp::Temperature,
    cpl,
    pms::AbstractVector,
    oa::AbstractVector,
    lambda_r1_r2_b::Tuple{Length,Length,Length},
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    digits::Integer,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    size::NTuple{2,Int}=(800, 600),
)
    f = Figure(; size)
    cr_info = "$(cr.metadata[:description]), critical phasematches for temperature: $(float(temp |> u"K")) ($(float(temp |> u"°C")))"
    lambda_info = "$(round(u"nm", lambda_r1_r2_b[1]; digits)) (λ_r1, $(hi_or_lo_r1_r2_b[1])) + $(round(u"nm", lambda_r1_r2_b[2]; digits)) (λ_r2, $(hi_or_lo_r1_r2_b[2])) = $(round(u"nm", lambda_r1_r2_b[3]; digits)) (λ_b, $(hi_or_lo_r1_r2_b[3]))"
    ax = Axis(
        f[1, 1],
        xlabel="ϕ",
        ylabel="θ",
        title=title = cr_info * "\n" * lambda_info
    )

    heatmap!(ax, ϕ_range .|> u"°", θ_range .|> u"°", all_delta_k'; colormap=COLORMAP_HEATMAP, colorrange=(-scale_limit, scale_limit), inspectable=false)

    plot_polar_contours!(ax, cpl, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_polar_pms!(ax, pms)

    if show_coordinates
        plot_polar_coordinate_markers!(ax)
    end

    if show_optical_axes
        plot_polar_optical_axes!(ax, oa)
    end

    # Legend(f[1, 2], ax)
    Colorbar(f[1, 2]; colormap=COLORMAP_HEATMAP, limits=(-scale_limit, scale_limit), label="Δk in µm⁻¹")
    DataInspector(ax)

    return f
end

function plot_pm_label(pm::CollinearPhaseMatch)
    buf = IOBuffer()
    print(buf, pm)
    return String(take!(buf))
end

function polar_plot_pm_label(
    plot,
    index,
    position,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Length,Length,Length}
)
    θ_pm = position[2] * u"°"
    ϕ_pm = position[1] * u"°"
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
    return plot_pm_label(pm)
end

function plot_polar_contours!(
    ax::Axis,
    cpl,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Length,Length,Length};
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θi_cl = [v[1] for v in cl.vertices] .|> u"°"
        ϕi_cl = [v[2] for v in cl.vertices] .|> u"°"

        label_fun = (plot, index, position) -> polar_plot_pm_label(plot, index, position, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)

        lines!(
            ax, ϕi_cl, θi_cl;
            linewidth=2, fxaa=true, inspectable=true,
            inspector_label=label_fun,
            label=(i == 1 ? "Phase matches" : nothing), color=COL_CONTOUR
        )
    end
end

function plot_polar_pms!(ax::Axis, pms::AbstractVector)
    for (i, pm) in enumerate(pms)
        scatter!(
            ax, pm.phi_pm |> u"°", pm.theta_pm |> u"°";
            markersize=10,
            label="Probed phase match $(i)", inspectable=true,
            inspector_label=(plot, index, position) -> ("Probed phase match $(i)\n" * plot_pm_label(pm))
        )
    end
end

function plot_polar_coordinate_markers!(ax::Axis)
    digits = 3
    all_coordinate_labels = [
        (0.0, 0.0, "+Z", (:left, :bottom)),
        (0.0, π / 2, "+Z", (:left, :bottom)),
        (0.0, π, "+Z", (:left, :bottom)),
        (0.0, 3π / 2, "+Z", (:left, :bottom)),
        (0.0, 2π, "+Z", (:right, :bottom)),
        (π / 2, 0.0, "+X", (:left, :bottom)),
        (π / 2, π / 2, "+Y", (:left, :bottom)),
        (π / 2, π, "-X", (:left, :bottom)),
        (π / 2, 3π / 2, "-Y", (:left, :bottom)),
        (π / 2, 2π, "+X", (:right, :bottom)),
        (π, 0.0, "-Z", (:left, :top)),
        (π, π / 2, "-Z", (:right, :top)),
        (π, π, "-Z", (:right, :top)),
        (π, 3π / 2, "-Z", (:right, :top)),
        (π, 2π, "-Z", (:right, :top)),
    ]
    for (θ, ϕ, lab, a) in all_coordinate_labels
        θ = θ |> u"°"
        ϕ = ϕ |> u"°"
        coord_label = "$(lab) axis\nθ = $(round(u"°", θ; digits)) ($(round(u"rad", θ; digits)))\nϕ = $(round(u"°", ϕ; digits)) ($(round(u"rad", ϕ; digits)))"
        scatter!(ax, ϕ, θ; color=COL_COORDS, inspectable=true, inspector_label=(plot, index, position) -> coord_label)
        text!(ax, ϕ, θ; text=lab, color=COL_COORDS, align=a)
    end
end

function plot_polar_optical_axes!(ax::Axis, oa)
    digits = 3
    colors = [COL_R1, COL_R2, COL_B]
    labels = ["λ_r1", "λ_r2", "λ_b"]

    for (angle, color, lab) in zip(oa, colors, labels)
        axis_label = "Optical axis $(lab)\n|θ| = $(round(u"°", angle; digits)) ($(round(u"rad", angle; digits)))\nϕ = $(0.0u"°") ($(0.0u"rad"))"
        scatter!(
            ax,
            [0, 0, π, π, 2π, 2π] .|> u"°",
            [angle, π - angle, angle, π - angle, angle, π - angle] .|> u"°";
            color, markersize=8, label=lab,
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
    temp::Temperature,
    cpl,
    pms::AbstractVector,
    oa::AbstractVector,
    lambda_r1_r2_b::Tuple{Length,Length,Length},
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    digits::Integer,
    axis_length::Real,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    size::NTuple{2,Int}=(800, 600),
)
    f = Figure(; size)
    cr_info = "$(cr.metadata[:description]), critical phasematches for temperature: $(float(temp |> u"K")) ($(float(temp |> u"°C")))"
    lambda_info = "$(round(u"nm", lambda_r1_r2_b[1]; digits)) (λ_r1, $(hi_or_lo_r1_r2_b[1])) + $(round(u"nm", lambda_r1_r2_b[2]; digits)) (λ_r2, $(hi_or_lo_r1_r2_b[2])) = $(round(u"nm", lambda_r1_r2_b[3]; digits)) (λ_b, $(hi_or_lo_r1_r2_b[3]))"
    ax = Axis3(
        f[1, 1],
        azimuth=0.1π, elevation=0.05π,
        aspect=:data, viewmode=:fit,
        title=cr_info * "\n" * lambda_info
    )
    hidedecorations!(ax, grid=false)

    plot_sphere_mesh!(ax, θ_range, ϕ_range, all_delta_k, scale_limit)
    plot_sphere_contours!(ax, cpl, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_sphere_pms!(ax, pms, axis_length)

    if show_coordinates
        plot_sphere_coordinates!(ax, axis_length)
    end

    if show_optical_axes
        plot_sphere_optical_axes!(ax, oa, axis_length)
    end

    # Legend(f[1, 2], ax)
    Colorbar(
        f[1, 2]; colormap=COLORMAP_HEATMAP,
        limits=(-scale_limit, scale_limit),
        label="Δk in µm⁻¹"
    )
    DataInspector(ax)

    return f
end

function plot_sphere_mesh!(ax::Axis3, θ_range, ϕ_range, all_delta_k, scale_limit)
    points = [Point3f(angles_to_vector(θ, ϕ)...) for θ in θ_range, ϕ in ϕ_range]
    faces = decompose(QuadFace{GLIndex}, Tessellation(Rect(0, 0, 1, 1), size(all_delta_k)))
    normals = normalize.(vec(points))
    gb_mesh = GeometryBasics.Mesh(vec(points), faces; normal=normals)
    mesh!(ax, gb_mesh, color=all_delta_k[:]; colormap=COLORMAP_HEATMAP, colorrange=(-scale_limit, scale_limit), inspectable=false)
end

function sphere_plot_pm_label(
    plot,
    index,
    position,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Length,Length,Length},
)
    θ_pm, ϕ_pm = vector_to_angles(position)
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
    return plot_pm_label(pm)
end

function plot_sphere_contours!(
    ax::Axis3,
    cpl,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    lambda_r1_r2_b::Tuple{Length,Length,Length}
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θ = [v[1] for v in cl.vertices] .|> u"°"
        ϕ = [v[2] for v in cl.vertices] .|> u"°"

        label_fun = (plot, index, position) -> sphere_plot_pm_label(plot, index, position, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)

        points = [Point3f(angles_to_vector(θi, ϕi)...) for (θi, ϕi) in zip(θ, ϕ)]
        lines!(ax, points; color=COL_CONTOUR, linewidth=4, label=(i == 1 ? "Phase matches" : nothing), inspectable=true, inspector_label=label_fun,
        )
    end
end

function plot_sphere_pms!(
    ax::Axis3,
    pms,
    axis_length
)
    fac_r1 = [1 + 2 * (axis_length - 1) / 3, axis_length]
    fac_r2 = [1 + (axis_length - 1) / 3, 1 + 2 * (axis_length - 1) / 3]
    fac_b = [1, 1 + (axis_length - 1) / 3]

    len_E_dir = (axis_length - 1) / 3

    lambda_labels = ["λ_r1", "λ_r2", "λ_b"]
    factors = [fac_r1, fac_r2, fac_b]
    colors = [COL_R1, COL_R2, COL_B]

    for (i, pm) in enumerate(pms)
        θi_pm = pm.theta_pm
        ϕi_pm = pm.phi_pm
        E_dirs = pm.E_dir_r1_r2_b

        ri_pm = Point3f(angles_to_vector(θi_pm, ϕi_pm)...)

        scatter!(ax, axis_length .* ri_pm; markersize=15, label="Probed phase match $(i)", inspectable=true,
            inspector_label=(plot, index, position) -> ("Probed phase match $(i)\n" * plot_pm_label(pm)))

        for (i, (lab, fac, E_dir, color)) in enumerate(zip(lambda_labels, factors, E_dirs, colors))
            label_fun = (plot, index, position) -> ("$(lab)")
            lines!(
                ax,
                [fac[1] .* ri_pm, fac[2] .* ri_pm];
                color, linewidth=4, inspectable=true, inspector_label=label_fun
            )

            scale = sum(fac) / 2
            lines!(
                ax,
                scale .* ri_pm[1] .+ len_E_dir .* E_dir[1] .* [-1, +1],
                scale .* ri_pm[2] .+ len_E_dir .* E_dir[2] .* [-1, +1],
                scale .* ri_pm[3] .+ len_E_dir .* E_dir[3] .* [-1, +1];
                color, linewidth=4, inspectable=true, inspector_label=label_fun
            )
        end
    end
end

function plot_sphere_coordinates!(ax::Axis3, axis_length)
    digits = 3
    coords = [
        ([-1, 1], [0, 0], [0, 0], ["-X", "+X"]),
        ([0, 0], [-1, 1], [0, 0], ["-Y", "+Y"]),
        ([0, 0], [0, 0], [-1, 1], ["-Z", "+Z"])
    ]
    for (x, y, z, label) in coords
        θ, ϕ = vector_to_angles([x[2], y[2], z[2]])
        θ = θ |> u"°"
        ϕ = ϕ |> u"°"
        coord_label = "$(label) axis\nθ = $(round(u"°", θ; digits)) ($(round(u"rad", θ; digits)))\nϕ = $(round(u"°", ϕ; digits)) ($(round(u"rad", ϕ; digits)))"
        lines!(ax, axis_length .* x, axis_length .* y, axis_length .* z; color=COL_COORDS, inspectable=true, inspector_label=(plot, index, position) -> coord_label, linewidth=3)
        text!(ax, axis_length .* x[2], axis_length .* y[2], axis_length .* z[2]; text=label[2], color=COL_COORDS)
        text!(ax, axis_length .* x[1], axis_length .* y[1], axis_length .* z[1]; text=label[1], color=COL_COORDS)
    end
end

function plot_sphere_optical_axes!(ax::Axis3, oa, axis_length)
    digits = 3
    colors = [COL_R1, COL_R2, COL_B]
    labels = ["λ_r1", "λ_r2", "λ_b"]

    for (angle, color, lab) in zip(oa, colors, labels)
        axis_label = "Optical axis $(lab)\n|θ| = $(round(u"°", angle; digits)) ($(round(u"rad", angle; digits)))\nϕ = $(0.0u"°") ($(0.0u"rad"))"
        dir_pos = axis_length * Point3f(sin(angle), 0, cos(angle))
        dir_neg = -dir_pos
        lines!(ax, [dir_neg, dir_pos]; color, linewidth=4, label="Optical axes " * lab, inspectable=true, inspector_label=(plot, index, position) -> axis_label)

        dir_pos_neg = axis_length * Point3f(sin(-angle), 0, cos(-angle))
        dir_neg_neg = -dir_pos_neg
        lines!(ax, [dir_neg_neg, dir_pos_neg]; color, linewidth=4, inspectable=true, inspector_label=(plot, index, position) -> axis_label)
    end
end