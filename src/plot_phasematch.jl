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

            # for i in eachindex(segments_cb)
            #     for p in eachindex(segments_cb[i])
            #         pm = find_nearest_pm_along_lambda_r_b(
            #             hi_or_lo_r1_r2_b,
            #             cr;
            #             lambda_r1=segments_cr[i][p],
            #             lambda_r2=nothing,
            #             lambda_b=segments_cb[i][p],
            #             temp,
            #             principal_axis=principal_axis,
            #         )
            #         if isnothing(pm)
            #             pop!(segments_cb[i][p])
            #             pop!(segments_cr[i][p])
            #         else
            #             segments_cb[i][p] = pm.lambda_r1_r2_b[3]
            #             segments_cr[i][p] = pm.lambda_r1_r2_b[1]
            #         end
            #     end
            # end

            # for i in eachindex(segments_cb)
            #     lines!(segments_cb[i], segments_cr[i])
            # end
            # !isempty(cb_intersections) && scatter!(cb_intersections, cr_intersections)

            # if !isempty(sign_switches)
            #     scatter!(cont_b[sign_switches], cont_r[sign_switches] .- 2 * cont_b[sign_switches])
            # end
            # lines!(cont_b, 1 ./ (1 ./ cont_b - 1 ./ cont_r))
        end
    end
    all_segments_b = reduce(vcat, all_segments_b)
    all_segments_r = reduce(vcat, all_segments_r)
    all_cb_intersections = reduce(vcat, all_cb_intersections)
    all_cr_intersections = reduce(vcat, all_cr_intersections)
    return all_segments_b, all_segments_r, all_cb_intersections, all_cr_intersections
end

# function find_sorted_segments(sorted_indices)
#     segments = []
#     visited = falses(length(sorted_indices))

#     for i in eachindex(sorted_indices)
#         if !visited[i]
#             segment = [i]
#             visited[i] = true

#             # Forward traversal
#             current = i
#             while true
#                 next_idx = findfirst(n -> !visited[n], sorted_indices[current])
#                 isnothing(next_idx) && break
#                 current = sorted_indices[current][next_idx]
#                 push!(segment, current)
#                 visited[current] = true
#             end

#             # Backward traversal
#             current = i
#             while true
#                 prev_point = nothing
#                 for neighbor in sorted_indices[current]
#                     if !visited[neighbor] && (current in sorted_indices[neighbor])
#                         prev_point = neighbor
#                         break
#                     end
#                 end
#                 isnothing(prev_point) && break
#                 prepend!(segment, prev_point)
#                 visited[prev_point] = true
#                 current = prev_point
#             end

#             push!(segments, segment)
#         end
#     end

#     return segments
# end

# function calc_raw_noncritical_pm_points(
#     principal_axis::Symbol,
#     hi_or_lo_r1_r2_b::AbstractVector{Symbol},
#     cr::NonlinearCrystal;
#     lambda_b_min::Union{Nothing,Length}=nothing,
#     lambda_b_max::Union{Nothing,Length}=nothing,
#     lambda_r12_min::Union{Nothing,Length}=nothing,
#     lambda_r12_max::Union{Nothing,Length}=nothing,
#     temp::Temperature=default_temp(cr),
#     ngrid=50,
#     tol=1e-14u"nm^-1",
# )
#     @assert principal_axis in [:X, :Y, :Z]
#     @assert ngrid >= 2 "The grid must contain at least two points."

#     isnothing(lambda_b_min) && (lambda_b_min = valid_lambda_range(cr)[1])
#     isnothing(lambda_b_max) && (lambda_b_max = valid_lambda_range(cr)[2] / 2)
#     isnothing(lambda_r12_min) && (lambda_r12_min = valid_lambda_range(cr)[1])
#     isnothing(lambda_r12_max) && (lambda_r12_max = valid_lambda_range(cr)[2])

#     range_lambda_b = LinRange(lambda_b_min, lambda_b_max, ngrid)
#     all_pms_b = typeof(lambda_b_min)[]
#     all_pms_r1 = typeof(lambda_b_min)[]
#     all_pms_r2 = typeof(lambda_b_min)[]
#     for lambda_b in range_lambda_b
#         pms = find_all_pms_along_dimension(hi_or_lo_r1_r2_b, cr; lambda_b_fixed=lambda_b, temp_min=temp, temp_max=temp, principal_axis, ngrid, tol)
#         for pm in pms
#             # if true#pm.lambda_r1_r2_b[2] < pm.lambda_r1_r2_b[1]
#             push!(all_pms_b, lambda_b)
#             push!(all_pms_r1, pm.lambda_r1_r2_b[1])
#             push!(all_pms_r2, pm.lambda_r1_r2_b[2])
#             # end
#         end
#     end

#     range_lambda_r1 = lambda_r12_min:step(range_lambda_b):lambda_r12_max
#     for lambda_r1 in range_lambda_r1
#         pms = find_all_pms_along_dimension(hi_or_lo_r1_r2_b, cr; lambda_r1_fixed=lambda_r1, temp_min=temp, temp_max=temp, principal_axis, ngrid, tol)
#         for pm in pms
#             # if true#pm.lambda_r1_r2_b[2] < lambda_r1
#             push!(all_pms_b, pm.lambda_r1_r2_b[3])
#             push!(all_pms_r1, lambda_r1)
#             push!(all_pms_r2, pm.lambda_r1_r2_b[2])
#             # end
#         end
#     end

#     # @show step(range_lambda_b)
#     # @show indices_within_d, dists_within_d = find_neighbors_within_distance(all_pms_b, all_pms_r1, (1 + 1e-4) * sqrt(2) * step(range_lambda_b)) 
#     # segments = find_sorted_segments(indices_within_d)

#     # all_dir_r1_b = [reverse(delta_k_gradient_r1_b(principal_axis, hi_or_lo_r1_r2_b, cr; lambda_r1=all_pms_r1[i], lambda_b=all_pms_b[i])) .* [1, -1] for i in eachindex(all_pms_b)]

#     # idx = sortperm(all_pms_b)
#     # all_pms_b = all_pms_b[idx]
#     # all_pms_r1 = all_pms_r1[idx]
#     # all_pms_r2 = all_pms_r2[idx]
#     return all_pms_b, all_pms_r1, all_pms_r2
# end

function plot_single_noncritical_pm(
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

    f = Figure(size=(1000, 800))
    uc = Makie.UnitfulConversion(u"µm"; units_in_label=true)
    ax = Axis(
        f[1, 1],
        xlabel="λ_b",
        ylabel="λ_r12",
        title="$(cr.metadata[:description])\nNoncritical phasematches along positive $(principal_axis) axis, polarization directions: $(hi_or_lo_r1_r2_b)", # TODO: Add lambda names and o/e
        dim1_conversion=uc,
        dim2_conversion=uc,
    )

    vline_lab = (plot, idx, pos) -> begin
        pos = pos * u"µm"
        marker_x[] = [ustrip(u"µm", pos[1])]
        marker_x_b[] = [pos[1]]
        is_r1 = pos[2] > 2 * pos[1]
        marker_y_r1[] = is_r1 ? [pos[2]] : [1 / (1 / pos[1] - 1 / pos[2])]
        marker_y_r2[] = is_r1 ? [1 / (1 / pos[1] - 1 / pos[2])] : [pos[2]]

        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_r1_r2_b, cr; lambda_r1=marker_y_r1[][1], lambda_b=marker_x_b[][1], temp, principal_axis, ngrid, tol)
        return plot_pm_label(pm)
    end

    shg_lab = (plot, idx, pos) -> begin
        pos = pos * u"µm"
        pm = find_nearest_pm_along_lambda_r_b(hi_or_lo_r1_r2_b, cr; lambda_r1=pos[2], lambda_b=pos[1], temp, principal_axis, ngrid, tol)
        return "SHG point\n" * plot_pm_label(pm)
    end

    vline_clear = (inspector, plot) -> begin
        marker_x[] = [NaN]
        marker_x_b[] = [NaN * u"µm"]
        marker_y_r1[] = [NaN * u"µm"]
        marker_y_r2[] = [NaN * u"µm"]
        return nothing
    end

    [lines!(ax, all_segments_b[i], all_segments_r1[i], color=COL_R1, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
    [lines!(ax, all_segments_b[i], 1 ./ (1 ./ all_segments_b[i] .- 1 ./ all_segments_r1[i]), color=COL_R2, linewidth=2, inspectable=true, inspector_label=vline_lab, inspector_clear=vline_clear) for i in eachindex(all_segments_b)]
    scatter!(ax, all_cb_intersections, all_cr_intersections; inspectable=true, inspector_label=shg_lab)

    marker_x = Observable([NaN]) # Workaround: vlines! is not yet compatible with unit observables
    marker_x_b = Observable([NaN * u"µm"])
    marker_y_r1 = Observable([NaN * u"µm"])
    marker_y_r2 = Observable([NaN * u"µm"])
    vlines!(ax, marker_x; inspectable=false, color=:gray, linewidth=2, alpha=0.5)
    scatter!(ax, marker_x_b, marker_y_r1; inspectable=false, color=COL_R1)
    scatter!(ax, marker_x_b, marker_y_r2; inspectable=false, color=COL_R2)

    DataInspector(f)

    return f
end

# function on_hover(inspector)
#     parent = inspector.root
#     (inspector.attributes.enabled[] && is_mouseinside(parent)) || return Consume(false)

#     mp = mouseposition_px(parent)
#     should_clear = true
#     for (plt, idx) in pick_sorted(parent, mp, inspector.attributes.range[])
#         if to_value(get(plt.attributes, :inspectable, true))
#             # show_data should return true if it created a tooltip
#             if show_data_recursion(inspector, plt, idx)
#                 should_clear = false
#                 break
#             end
#         end
#     end

#     if should_clear
#         plot = inspector.selection
#         if to_value(get(plot, :inspector_clear, automatic)) !== automatic
#             plot[:inspector_clear][](inspector, plot)
#         end
#         inspector.plot.visible[] = false
#         inspector.attributes.indicator_visible[] = false
#         inspector.plot.offset.val = inspector.attributes.offset[]
#     end

#     return Consume(false)
# end

function plot_single_noncritial_pm!(
    ax::Axis,
    θ::Angle,
    ϕ::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal,
    range_lambda_b,
    range_lambda_r12,
    temp::Temperature=default_temp(cr),
)
    col_r12 = (colormap=:vik10, colorrange=(1, 10))
    color_r1 = 9
    color_r2 = 8

    # Only search along positive direction of the specified axes
    all_delta_k = zeros(length(range_lambda_b), length(range_lambda_r12))
    for i in CartesianIndices(all_delta_k)
        (i_b, i_r12) = (i[1], i[2])
        lambda_b = range_lambda_b[i_b]
        lambda_r12 = range_lambda_r12[i_r12]
        lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_b, lambda_r1=lambda_r12)
        if is_lambda_valid(lambda_r2, cr)
            all_delta_k[i_b, i_r12] = ustrip(u"m^-1", delta_k(θ, ϕ, hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b))
        else
            all_delta_k[i_b, i_r12] = NaN
        end
    end

    # scale_limit = maximum(abs.(all_delta_k[.!isnan.(all_delta_k)]))
    # heatmap!(ax, range_lambda_b, range_lambda_r12, all_delta_k; colormap=:vik10, colorrange=(-scale_limit, scale_limit))

    # Extract Δk=0 isolines using Makie contour function
    cpl = Makie.Contours.contour(range_lambda_b, range_lambda_r12, all_delta_k, 0.0)
    for cl in cpl.lines
        cont_b = ([v[1] for v in cl.vertices])
        cont_r = ([v[2] for v in cl.vertices])
        lines!(ax, cont_b, cont_r; color=color_r1, col_r12..., linewidth=4, label="λ_r1 for all noncritial phasematches", inspectable=true)
    end
    # contour!(ax, range_lambda_b, range_lambda_r12, all_delta_k, levels=[0.0])

    # data_valid = .!isnan.(cont_r) .&& .!isnan.(cont_b)
    # (isempty(cont_r[data_valid]) || isempty(cont_b[data_valid])) && continue
    # cont_b = cont_b[data_valid]
    # cont_r = cont_r[data_valid]

    # Split cont_r lines in upper red line (r1) and lower red line (r2)
    # cont_b_with_shg = []
    # cont_r_with_shg = []
    # signs_diff_r12 = []
    # shg_indices = []
    # last_diff_r12 = nothing
    # for i in eachindex(cont_b)
    #     diff_r12 = 2 * cont_b[i] - cont_r[i]

    #     # Add SHG point when difference between red curves switches signs
    #     if !isnothing(last_diff_r12) && sign(diff_r12) != sign(last_diff_r12)
    #         frac = -last_diff_r12 / (diff_r12 - last_diff_r12)
    #         cont_b_shg = cont_b[i-1] + (cont_b[i] - cont_b[i-1]) * frac
    #         cont_r12_shg = cont_r[i-1] + (cont_r[i] - cont_r[i-1]) * (cont_b_shg - cont_b[i-1]) / (cont_b[i] - cont_b[i-1])
    #         push!(signs_diff_r12, 0)
    #         push!(cont_r_with_shg, cont_r12_shg)
    #         push!(cont_b_with_shg, cont_b_shg)
    #         push!(shg_indices, length(cont_b_with_shg))
    #     end

    #     # Add usual point and note sign as an indicator if the curve is r1 or r2
    #     push!(signs_diff_r12, sign(diff_r12))
    #     push!(cont_r_with_shg, cont_r[i])
    #     push!(cont_b_with_shg, cont_b[i])
    #     last_diff_r12 = diff_r12
    # end

    # # Plot all line segments for r1 and r2 separated by shg_points
    # label_fun = (plot, index, position; lambda_r_switch=:lambda_r1) -> noncritical_pm_label(
    #     plot,
    #     index,
    #     position,
    #     θ,
    #     ϕ,
    #     cr,
    #     temp,
    #     hi_or_lo_r1_r2_b,
    #     lambda_r_switch,
    # )

    # cont_r_other_with_shg = (1 ./ (1 ./ cont_b_with_shg .- 1 ./ cont_r_with_shg))
    # # i_segs = [1; shg_indices; length(cont_b_with_shg)]
    # # signs_diff_r12[i_segs] # TODO!

    # last_s = 1
    # for s in [shg_indices; length(cont_b_with_shg)]
    #     if signs_diff_r12[last_s+1] <= 0
    #         lines!(ax, cont_b_with_shg[last_s:s] .* u"µm", cont_r_with_shg[last_s:s] .* u"µm"; color=color_r1, col_r12..., linewidth=4, label="λ_r1 for all noncritial phasematches", inspectable=true, inspector_label=(args...) -> label_fun(args...; lambda_r_switch=:lambda_r1))
    #         lines!(ax, cont_b_with_shg[last_s:s] .* u"µm", cont_r_other_with_shg[last_s:s] .* u"µm"; color=color_r2, col_r12..., linewidth=4, label="λ_r2 for all noncritial phasematches", inspectable=true, inspector_label=(args...) -> label_fun(args...; lambda_r_switch=:lambda_r2))
    #     else
    #         lines!(ax, cont_b_with_shg[last_s:s] .* u"µm", cont_r_with_shg[last_s:s] .* u"µm"; color=color_r2, col_r12..., linewidth=4, label="λ_r2 for all noncritial phasematches", inspectable=true, inspector_label=(args...) -> label_fun(args...; lambda_r_switch=:lambda_r2))
    #         lines!(ax, cont_b_with_shg[last_s:s] .* u"µm", cont_r_other_with_shg[last_s:s] .* u"µm"; color=color_r1, col_r12..., linewidth=4, label="λ_r1 for all noncritial phasematches", inspectable=true, inspector_label=(args...) -> label_fun(args...; lambda_r_switch=:lambda_r1))
    #     end

    #     last_s = s
    # end
    # scatter!(ax, cont_b_with_shg[shg_indices] .* u"µm", cont_r_with_shg[shg_indices] .* u"µm"; markersize=10, label="SHG noncritial phasematch", inspectable=true, inspector_label=label_fun)

    # lines!(ax, cont_b .* u"µm", cont_r .* u"µm"; color=color_r1, col_r12..., linewidth=4, label="λ_r1 for all noncritial phasematches", inspectable=true)
    # lines!(ax, cont_b .* u"µm", 1 ./ (1 ./ cont_b .- 1 ./ cont_r) .* u"µm"; color=color_r1, col_r12..., linewidth=4, label="λ_r1 for all noncritial phasematches", inspectable=true)
    # end
    return nothing
end

function plot_noncritical_pms(
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    all_temp::AbstractVector{<:Temperature}=[default_temp(cr)],
    ngrid=100,
    # tol=1e-14u"nm^-1",
)
    isnothing(lambda_b_min) && (lambda_b_min = valid_lambda_range(cr)[1])
    isnothing(lambda_b_max) && (lambda_b_max = valid_lambda_range(cr)[2] / 2)
    isnothing(lambda_r12_min) && (lambda_r12_min = valid_lambda_range(cr)[1])
    isnothing(lambda_r12_max) && (lambda_r12_max = valid_lambda_range(cr)[2])
    range_lambda_b = LinRange(lambda_b_min, lambda_b_max, ngrid) .|> u"µm"
    range_lambda_r12 = LinRange(lambda_r12_min, lambda_r12_max, ngrid) .|> u"µm"

    axes = [:X, :Y, :Z]

    f = Figure(size=(1500, 500))
    all_ax = []
    for (ia, a) in enumerate(axes)
        uc = Makie.UnitfulConversion(u"µm"; units_in_label=true)
        push!(
            all_ax,
            Axis(
                f[1, ia],
                xlabel="λ_b",
                ylabel="λ_r12",
                title="Along positive $(a) axis, Polarization direction: $(hi_or_lo_r1_r2_b)", # TODO: Add lambda names and o/e
                dim1_conversion=uc,
                dim2_conversion=uc,
            )
        )

        for temp in all_temp
            plot_single_noncritial_pm!(all_ax[end], axes_to_θ_ϕ(a)[1]..., hi_or_lo_r1_r2_b, cr, range_lambda_b, range_lambda_r12, temp)
        end

        DataInspector(all_ax[end])
    end

    linkxaxes!(all_ax...)
    linkyaxes!(all_ax...)
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
)
    # Color settings
    col_coordinates = (color=5, colormap=:Greys, colorrange=(1, 10))
    col_contour = (color=1, colormap=:Set1_9, colorrange=(1, 9))
    col_r1 = (color=9, colormap=:vik10, colorrange=(1, 10))
    col_r2 = (color=8, colormap=:vik10, colorrange=(1, 10))
    col_b = (color=3, colormap=:vik10, colorrange=(1, 10))
    col_heatmap = (colormap=:vik,)

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

    plot_polar_contours!(ax, cpl, col_contour, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_polar_pms!(ax, pms)

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
    col_contour,
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
            color=:black, linewidth=2, fxaa=true, inspectable=true,
            inspector_label=label_fun,
            label=(i == 1 ? "Phase matches" : nothing), col_contour...
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
    plot_sphere_contours!(ax, cpl, col_contour, cr, temp, hi_or_lo_r1_r2_b, lambda_r1_r2_b)
    plot_sphere_pms!(ax, pms, axis_length; col_r1, col_r2, col_b)

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
    col_contour,
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
        lines!(ax, points; col_contour..., linewidth=4, label=(i == 1 ? "Phase matches" : nothing), inspectable=true, inspector_label=label_fun,
        )
    end
end

function plot_sphere_pms!(
    ax::Axis3,
    pms,
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

    for (i, pm) in enumerate(pms)
        θi_pm = pm.theta_pm
        ϕi_pm = pm.phi_pm
        E_dirs = pm.E_dir_r1_r2_b

        ri_pm = Point3f(angles_to_vector(θi_pm, ϕi_pm)...)

        scatter!(ax, axis_length .* ri_pm; markersize=15, label="Probed phase match $(i)", inspectable=true,
            inspector_label=(plot, index, position) -> ("Probed phase match $(i)\n" * plot_pm_label(pm)))

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