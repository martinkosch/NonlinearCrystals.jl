export calc_noncritical_pm_lines, plot_critical_pms, plot_noncritical_pms, plot_delta_k_map

function noncritical_pm_label(
    plot,
    index,
    position,
    θ_pm,
    ϕ_pm,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_rrb::NTuple{3,Symbol},
    lambda_r_switch::Symbol,
)
    lambda_b = position[1] * u"µm"
    lambda_r1 = lambda_r2 = nothing
    if lambda_r_switch == :lambda_r1
        lambda_r1 = position[2] * u"µm"
    elseif lambda_r_switch == :lambda_r2
        lambda_r2 = position[2] * u"µm"
    end
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_rrb, cr; temp, lambda_r1, lambda_r2, lambda_b)
    return isnothing(pm) ? "" : plot_pm_label(pm)
end

function calc_raw_noncritical_pm_lines(
    principal_axis::Symbol,
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
)
    @assert principal_axis in [:X, :Y, :Z]
    θ_ϕ_principal_axis = axes_to_θ_ϕ(principal_axis)[1]

    range_lambda_r12, range_lambda_b, all_delta_k = calc_delta_k_map(
        θ_ϕ_principal_axis...,
        hi_or_lo_rrb,
        cr;
        lambda_b_min,
        lambda_b_max,
        lambda_r12_min,
        lambda_r12_max,
        temp,
        ngrid,
    )

    # Extract Δk=0 isolines using Makie contour function
    cpl = Makie.Contours.contour(
        range_lambda_b,
        range_lambda_r12,
        ustrip.(u"m^-1", all_delta_k),
        0.0,
        VT=NTuple{2,eltype(range_lambda_b)}
    )
    return cpl # Makie.jl contour object with list of `lines` containing list of `vertices`
end

function plot_delta_k_heatmap(
    principal_axis::Symbol,
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
    size::NTuple{2,Int}=(800, 600),
)
    @assert principal_axis in [:X, :Y, :Z]
    θ_ϕ_principal_axis = axes_to_θ_ϕ(principal_axis)[1]

    range_lambda_r12, range_lambda_b, all_delta_k = calc_delta_k_map(
        θ_ϕ_principal_axis...,
        hi_or_lo_rrb,
        cr;
        lambda_b_min,
        lambda_b_max,
        lambda_r12_min,
        lambda_r12_max,
        temp,
        ngrid,
    )

    f = Figure(; size)
    uc = Makie.UnitfulConversion(u"µm"; units_in_label=true)
    ax = Axis(
        f[1, 1],
        xlabel="λ_b",
        ylabel="λ_r12",
        title="$(cr.metadata[:description])\nΔk along positive $(principal_axis) axis, temperature: $(round(u"K", temp; digits=3)) ($(round(u"°C", temp; digits=3)))", # , polarization directions: $(hi_or_lo_rrb) # TODO: Add lambda names and o/e
        dim1_conversion=uc,
        dim2_conversion=uc,
    )
    scale_limit = ustrip(u"m^-1", maximum(abs.(all_delta_k[.!isnan.(all_delta_k)])))

    heatmap!(ax, range_lambda_b, range_lambda_r12, ustrip.(u"m^-1", all_delta_k); colorrange=(-scale_limit, scale_limit), colormap=:vik)
    contour!(ax, range_lambda_b, range_lambda_r12, ustrip.(u"m^-1", all_delta_k); levels=[0.0], color=COL_CONTOUR, linewidth=2)
    return f
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

"""
    calc_noncritical_pm_lines(principal_axis, hi_or_lo_rrb, cr; ...)

Computes the coordinates of Δk = 0 isolines in the λ_b-λ_r1 parameter space where **noncritical phasematching** is achieved
(i.e., phase matching occurs at fixed propagation direction, along the crystal's `principal_axis`).

Returns:
- `segments_b`: vector of λ_b curves (x-axis)
- `segments_r`: corresponding λ_r1 curves (y-axis)
- `cb_intersections`, `cr_intersections`: SHG symmetry intersection points (e.g., λ_b = λ_r / 2)
"""
function calc_noncritical_pm_lines(
    principal_axis::Symbol,
    hi_or_lo_rrb::NTuple{3,Symbol},
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
        hi_or_lo_rrb,
        cr;
        lambda_b_min,
        lambda_b_max,
        lambda_r12_min,
        lambda_r12_max,
        temp,
        ngrid
    )

    # If both red waves are polarized equally (type 1 phasematching), the plot lines show a helpful phase match symmetry around the line cont_r_raw = 2 * cont_b_raw
    is_type_one = (hi_or_lo_rrb[1] == hi_or_lo_rrb[2])

    all_segments_r = []
    all_segments_b = []
    all_cb_intersections = []
    all_cr_intersections = []
    # Iterate over all identified raw contour lines and postprocess/refine them
    for cl in cpl.lines
        cont_b_raw = [v[1] for v in cl.vertices]
        cont_r_raw = [v[2] for v in cl.vertices]

        # Split contours into segments at all NaNs and remove all NaNs
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
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=50,
    tol=1e-14u"nm^-1",
)
    is_type_1 = hi_or_lo_rrb[1] === hi_or_lo_rrb[2]

    all_segments_b, all_segments_r1, all_cb_intersections, all_cr_intersections = calc_noncritical_pm_lines(
        principal_axis,
        hi_or_lo_rrb,
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

        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_rrb, cr; lambda_r1=marker_y_r1[][1], lambda_b=pos[1], temp, principal_axis, ngrid, tol)
        return "$(is_type_1 ? "Type 1" : "Type 2") noncritical phasematch\n" * plot_pm_label(pm)
    end

    shg_lab = (plot, idx, pos) -> begin
        pos = pos * u"µm"
        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_rrb, cr; lambda_r1=pos[2], lambda_b=pos[1], temp, principal_axis, ngrid, tol)
        return "$(is_type_1 ? "Type 1" : "Type 2") noncritical phasematch SHG point\n" * plot_pm_label(pm)
    end

    vline_clear = (inspector, plot) -> begin
        marker_x[] = [NaN]
        marker_x_b[] = [NaN * u"µm"]
        marker_y_r1[] = [NaN * u"µm"]
        marker_y_r2[] = [NaN * u"µm"]
        return nothing
    end

    [lines!(ax, all_segments_b[i], all_segments_r1[i]; linestyle=hi_or_lo_rrb[1] === :hi ? :dash : :solid, color=COL_R1, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
    [lines!(ax, all_segments_b[i], 1 ./ (1 ./ all_segments_b[i] .- 1 ./ all_segments_r1[i]); linestyle=hi_or_lo_rrb[2] === :hi ? :dash : :solid, color=COL_R2, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
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

"""
    plot_noncritical_pms(principal_axis, cr; ...)

Generates an interactive plot of all **noncritical phasematching lines** (Δk = 0) in λ_b-λ_r1 space for a given crystal `cr`
along the selected `principal_axis`.

Each plotted contour represents a set of wavelength triplets (λ_r1, λ_r2, λ_b) that satisfy phasematching without angular adjustment.
Supports both Type I (equal polarization) and Type II (cross-polarized) combinations.

You can control the polarization combination via `hi_or_lo_rrb`, limit the scan ranges, and adjust the resolution with `ngrid`.

Returns a GLMakie `Figure`.
"""
function plot_noncritical_pms(
    principal_axis::Symbol,
    cr::NonlinearCrystal;
    hi_or_lo_rrb::Union{NTuple{3,Symbol},AbstractVector{<:NTuple{3,Symbol}},Nothing}=nothing,
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
        title="$(cr.metadata[:description])\nNoncritical phasematches along positive $(principal_axis) axis, temperature: $(round(u"K", temp; digits=3)) ($(round(u"°C", temp; digits=3)))", # , polarization directions: $(hi_or_lo_rrb) # TODO: Add lambda names and o/e
        dim1_conversion=uc,
        dim2_conversion=uc,
    )

    if isnothing(hi_or_lo_rrb)
        hi_or_lo_rrb = bool_permutations(:hi, :lo, 3)
        setdiff!(hi_or_lo_rrb, [(:hi, :lo, :lo), (:hi, :lo, :hi), (:lo, :lo, :lo), (:hi, :hi, :hi)]) # Prevent double plotting of Type 2 noncritical phasematches and (usually) unphysical phasematch combinations
    elseif typeof(hi_or_lo_rrb) <: NTuple{3,Symbol}
        hi_or_lo_rrb = [hi_or_lo_rrb]
    end

    for hl in hi_or_lo_rrb
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

"""
    plot_critical_pms(cr::NonlinearCrystal; kwargs...) -> Figure

Visualizes and compares **critical phase-matching solutions** for a given nonlinear crystal `cr` across all possible
polarization configurations and propagation directions. This plot gives insight into how the phasematching 
characteristics vary with angle for a **fixed wavelength triplet** (λ_r1, λ_r2, λ_b).

Each horizontal segment in the figure corresponds to a distinct type of phase-matching configuration, and
the plot shows how key quantities vary along that solution contour (i.e., as a function of θ and ϕ for which Δk = 0).

#### Plotted Quantities per Phasematch
For each matched solution (Δk = 0), the following attributes are visualized:
- **Phase velocity / c₀** — the normalized refractive indices (n)
- **Group velocity / c₀** — inverse group indices (group velocity dispersion)
- **GDD** — group delay dispersion (β₂), unit: fs²/mm
- **Walkoff angle** — spatial beam walkoff, unit: mrad
- **ω BW × L** — angular frequency bandwidth product (Δω·L), from group velocity mismatch
- **T BW × L** — temperature tolerance (ΔT·L), from thermal dispersion mismatch
- **ϕ BW × L**, **θ BW × L** — angular acceptance bandwidths
- **|d_eff|** — effective nonlinearity (with Miller scaling applied)
- **ϕ**, **θ** — propagation angles
- **Type** — phase-matching type and polarization roles (label only)

#### Optional Arguments
- `hi_or_lo_rrb`: One or more polarization configurations (e.g., `(:hi, :hi, :lo)`)
- `lambda_r1`, `lambda_r2`, `lambda_b`: At least two must be specified, the third is inferred
- `temp`: Temperature at which to evaluate phasematching
- `n_points`: Angular resolution
- `size`: Plot size as a tuple (width, height)

#### Output
Returns a vertically stacked GLMakie `Figure` containing multiple linked plots, each showing one of the quantities above.
Each horizontal span corresponds to a continuous critical phase-matching curve (i.e., varying θ and ϕ for fixed λ and T).
"""
function plot_critical_pms(cr::NonlinearCrystal;
    hi_or_lo_rrb::Union{NTuple{3,Symbol},AbstractVector{<:NTuple{3,Symbol}},Nothing}=nothing,
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    n_points::Integer=100,
    skip_symmetric_pms::Bool=false,
    size::NTuple{2,Int}=(700, 1200)
)
    lambda_rrb = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    # Set hi/lo permutations
    if isnothing(hi_or_lo_rrb)
        hi_or_lo_rrb = bool_permutations(:hi, :lo, 3)
        setdiff!(hi_or_lo_rrb, [(:lo, :lo, :lo), (:hi, :hi, :hi)]) # Remove (usually) unphysical phasematch combinations
    elseif isa(hi_or_lo_rrb, NTuple{3,Symbol})
        hi_or_lo_rrb = [hi_or_lo_rrb]
    end

    # Compute phase matches
    all_pm = map(hi_or_lo_rrb) do hilo
        θ_range, ϕ_range, all_delta_k = compute_delta_k_grid(
            cr,
            hilo,
            lambda_rrb...,
            temp,
            n_points;
            theta_range_max=(skip_symmetric_pms ? 90u"°" : 180u"°"),
            phi_range_max=(skip_symmetric_pms ? 90u"°" : 360u"°"),
        )

        # Generate contours
        cpl = Makie.Contours.contours(ustrip.(u"rad", θ_range), ustrip.(u"rad", ϕ_range), all_delta_k, [0.0])

        return map(cpl.contours[1].lines) do cl
            return map(cl.vertices) do v
                return find_nearest_pm_along_theta_phi(
                    v[1] |> u"°",
                    v[2] |> u"°",
                    hilo,
                    cr;
                    temp,
                    lambda_r1=lambda_rrb[1],
                    lambda_r2=lambda_rrb[2],
                    lambda_b=lambda_rrb[3]
                )
            end
        end
    end

    # Set up data, units, and data sources
    axis_data = [
        ("Phase vel. / c₀", nothing, nothing, pm -> pm.n_rrb),
        ("Group vel. / c₀", nothing, nothing, pm -> pm.group_index_rrb),
        ("GDD", u"fs^2/mm", nothing, pm -> pm.beta2_rrb),
        ("Walkoff angle", u"mrad", nothing, pm -> pm.walkoff_angle_rrb),
        ("ω BW × L", u"GHz * cm", (0u"GHz * cm", 1000u"GHz * cm"), pm -> pm.bw_data.omega_L_bw),
        ("T BW × L", u"K * cm", (0u"K * cm", 100u"K * cm"), pm -> pm.bw_data.temp_L_bw),
        ("ϕ BW × L", u"mrad * cm", (0u"mrad * cm", 100u"mrad * cm"), pm -> pm.bw_data.phi_L_bw),
        ("θ BW × L", u"mrad * cm", (0u"mrad * cm", 100u"mrad * cm"), pm -> pm.bw_data.theta_L_bw),
        # ("|d_eff_no_miller|", u"pm/V", nothing, pm -> abs(pm.eff_data.d_eff_no_miller)),
        ("|d_eff|", u"pm/V", nothing, pm -> abs(pm.eff_data.d_eff)),
        ("ϕ", u"°", nothing, pm -> pm.phi_pm),
        ("θ", u"°", nothing, pm -> pm.theta_pm),
        ("Type", nothing, nothing, nothing),
    ]

    # Create figure and axes
    f = Figure(; size)
    axes = [
        Axis(f[i, 1],
            ylabel=ad[1],
            dim2_conversion=isnothing(ad[2]) ? nothing : Makie.UnitfulConversion(ad[2]; units_in_label=true))
        for (i, ad) in enumerate(axis_data)
    ]

    foreach(hidexdecorations!, axes[1:end-1])
    hidedecorations!(axes[end])
    linkxaxes!(axes)

    # Plot data
    start = 0
    for hl in all_pm, l in hl
        idx_range = start .+ (0:length(l)-1)

        for i in eachindex(axis_data)
            ax = axes[i]
            un = axis_data[i][2]
            lims = axis_data[i][3]
            extr = axis_data[i][4]
            isnothing(extr) && continue
            data = [extr(pm) for pm in l]
            if length(data[1]) == 3
                for (comp, col) in zip(1:3, [COL_R1, COL_R2, COL_B])
                    lines!(ax, idx_range, getindex.(data, comp), color=col, linewidth=2)
                end
            else
                lines!(ax, idx_range, data, linewidth=2)
            end
        end

        # Type labels
        vspan!(axes[end], [start], [start + length(l) - 1])

        label_text = if isa(cr, UnidirectionalCrystal)
            types = l[1].pm_data.pm_type[1]
            polars = ["$(l[1].hi_or_lo_rrb[i]) ($(types.o_or_e_rrb[i]))" for i in 1:3]
            "Type $(types.type) PM:\n$(join(polars, ", "))"
        else
            join(l[1].hi_or_lo_rrb, ", ")
        end
        text!(axes[end], start + (length(l) - 1) / 2, 0.5; text=label_text, align=(:center, :center))

        start += length(l) + 30
    end

    # Restrict automatic y axis limits to specified limits 
    for hl in all_pm, l in hl
        for i in eachindex(axis_data)
            ax = axes[i]
            un = axis_data[i][2]
            lims = axis_data[i][3]
            autolimits!(ax)
            if !isnothing(lims)
                current_ylims = ax.yaxis.attributes.limits[] .* un
                ylims!(
                    ax,
                    (
                        max(current_ylims[1], lims[1]),
                        min(current_ylims[2], lims[2])
                    ))
            end
        end
    end

    return f
end

"""
    plot_delta_k_map(hi_or_lo_rrb, cr; ...)

Plots a 2D or 3D map of the phasemismatch Δk over the angular domain (θ, ϕ) for a given polarization assignment (`hi_or_lo_rrb`)
and wavelength triplet.

Depending on `plot_type`:
- `:polar` → Shows Mercator style projected θ–ϕ map with Δk color-coded and phase match contours overlaid
- `:sphere` → Visualizes Δk on the surface of a 3D unit sphere using ray direction vectors
"""
function plot_delta_k_map(
    hi_or_lo_rrb::NTuple{3,Symbol},
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
    lambda_rrb = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    θ_range, ϕ_range, all_delta_k = compute_delta_k_grid(cr, hi_or_lo_rrb, lambda_rrb..., temp, n_points)
    scale_limit = maximum(abs.(all_delta_k))

    # Generate contours
    cpl = Makie.Contours.contours(ustrip.(u"rad", θ_range), ustrip.(u"rad", ϕ_range), all_delta_k, [0.0])

    # Compute optical axes
    oa = [optical_axis_angle(cr, λ, temp) for λ in lambda_rrb]

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
            lambda_rrb,
            hi_or_lo_rrb,
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
            lambda_rrb,
            hi_or_lo_rrb,
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
    return plot_delta_k_map(pm.hi_or_lo_rrb, pm.cr; lambda_r1=pm.lambda_rrb[1], lambda_r2=pm.lambda_rrb[2], lambda_b=pm.lambda_rrb[3], temp=pm.temp, pms=[pm], n_points, axis_length, digits, plot_type, show_coordinates, show_optical_axes)
end

function compute_delta_k_grid(
    cr::NonlinearCrystal,
    hi_or_lo_rrb,
    lambda_r1,
    lambda_r2,
    lambda_b,
    temp,
    n_points;
    theta_range_max::Angle=180u"°",
    phi_range_max::Angle=360u"°",
)
    θ_range = LinRange(0.0u"°", theta_range_max, n_points * cld(theta_range_max, 180u"°"))
    ϕ_range = LinRange(0.0u"°", phi_range_max, n_points * cld(phi_range_max, 180u"°"))

    # TODO: Use symmetry information for speedup
    if isa(cr, UnidirectionalCrystal)
        # No ϕ dependence, compute only one value and repeat
        all_delta_k = [
            delta_k(θ, 0.0u"°", hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp)
            for θ in θ_range
        ]
        all_delta_k = repeat(reshape(all_delta_k, :, 1), 1, length(ϕ_range))
    else
        all_delta_k = [delta_k(θ, ϕ, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in θ_range, ϕ in ϕ_range]
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
    lambda_rrb::Tuple{Length,Length,Length},
    hi_or_lo_rrb::NTuple{3,Symbol},
    digits::Integer,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    size::NTuple{2,Int}=(800, 600),
)
    f = Figure(; size)

    uc = Makie.UnitfulConversion(u"°"; units_in_label=true)
    cr_info = "$(cr.metadata[:description]), critical phasematches for temperature: $(round(u"K", temp; digits=3)) ($(round(u"°C", temp; digits=3)))"
    lambda_info = "$(round(u"nm", lambda_rrb[1]; digits)) (λ_r1, $(hi_or_lo_rrb[1])) + $(round(u"nm", lambda_rrb[2]; digits)) (λ_r2, $(hi_or_lo_rrb[2])) = $(round(u"nm", lambda_rrb[3]; digits)) (λ_b, $(hi_or_lo_rrb[3]))"
    ax = Axis(
        f[1, 1],
        xlabel="ϕ",
        ylabel="θ",
        title=title = cr_info * "\n" * lambda_info,
        dim1_conversion=uc,
        dim2_conversion=uc,
    )

    heatmap!(ax, ϕ_range .|> u"°", θ_range .|> u"°", all_delta_k'; colormap=COLORMAP_HEATMAP, colorrange=(-scale_limit, scale_limit), inspectable=false)

    plot_polar_contours!(ax, cpl, cr, temp, hi_or_lo_rrb, lambda_rrb)
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
    hi_or_lo_rrb::NTuple{3,Symbol},
    lambda_rrb::Tuple{Length,Length,Length}
)
    θ_pm = position[2] * u"°"
    ϕ_pm = position[1] * u"°"
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_rrb, cr; temp, lambda_r1=lambda_rrb[1], lambda_r2=lambda_rrb[2], lambda_b=lambda_rrb[3])
    return plot_pm_label(pm)
end

function plot_polar_contours!(
    ax::Axis,
    cpl,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_rrb::NTuple{3,Symbol},
    lambda_rrb::Tuple{Length,Length,Length};
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θi_cl = [v[1] for v in cl.vertices] .|> u"°"
        ϕi_cl = [v[2] for v in cl.vertices] .|> u"°"

        label_fun = (plot, index, position) -> polar_plot_pm_label(plot, index, position, cr, temp, hi_or_lo_rrb, lambda_rrb)

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
    lambda_rrb::Tuple{Length,Length,Length},
    hi_or_lo_rrb::NTuple{3,Symbol},
    digits::Integer,
    axis_length::Real,
    show_coordinates::Bool,
    show_optical_axes::Bool,
    size::NTuple{2,Int}=(800, 600),
)
    f = Figure(; size)
    cr_info = "$(cr.metadata[:description]), critical phasematches for temperature: $(round(u"K", temp; digits=3)) ($(round(u"°C", temp; digits=3)))"
    lambda_info = "$(round(u"nm", lambda_rrb[1]; digits)) (λ_r1, $(hi_or_lo_rrb[1])) + $(round(u"nm", lambda_rrb[2]; digits)) (λ_r2, $(hi_or_lo_rrb[2])) = $(round(u"nm", lambda_rrb[3]; digits)) (λ_b, $(hi_or_lo_rrb[3]))"
    ax = Axis3(
        f[1, 1],
        azimuth=0.1π, elevation=0.05π,
        aspect=:data, viewmode=:fit,
        title=cr_info * "\n" * lambda_info
    )
    hidedecorations!(ax, grid=false)

    plot_sphere_mesh!(ax, θ_range, ϕ_range, all_delta_k, scale_limit)
    plot_sphere_contours!(ax, cpl, cr, temp, hi_or_lo_rrb, lambda_rrb)
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
    hi_or_lo_rrb::NTuple{3,Symbol},
    lambda_rrb::Tuple{Length,Length,Length},
)
    θ_pm, ϕ_pm = vector_to_angles(position)
    pm = find_nearest_pm_along_theta_phi(θ_pm, ϕ_pm, hi_or_lo_rrb, cr; temp, lambda_r1=lambda_rrb[1], lambda_r2=lambda_rrb[2], lambda_b=lambda_rrb[3])
    return plot_pm_label(pm)
end

function plot_sphere_contours!(
    ax::Axis3,
    cpl,
    cr::NonlinearCrystal,
    temp::Temperature,
    hi_or_lo_rrb::NTuple{3,Symbol},
    lambda_rrb::Tuple{Length,Length,Length}
)
    for (i, cl) in enumerate(cpl.contours[1].lines)
        θ = [v[1] for v in cl.vertices] .|> u"°"
        ϕ = [v[2] for v in cl.vertices] .|> u"°"

        label_fun = (plot, index, position) -> sphere_plot_pm_label(plot, index, position, cr, temp, hi_or_lo_rrb, lambda_rrb)

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
        E_dirs = pm.E_dir_rrb

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