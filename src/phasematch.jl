export PhaseMatch, CollinearPhaseMatch, pm_wavelengths, find_all_pms_along_dimension, find_nearest_pm_along_lambda_r_b, find_nearest_pm_along_theta_phi, delta_k

abstract type PhaseMatch end

struct CollinearPhaseMatch{LT,TT,OT,AT,IT,WT,DT,ET,ST,RT,GT,B2T,B3T,FT,FS,LBW,TBW,HBW,PBW,CT} <: PhaseMatch
    lambda_r1_r2_b::Vector{<:LT}
    temp::TT
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
    S_0::FS
    lambda_L_bw::LBW
    temp_L_bw::TBW
    theta_L_bw::HBW
    phi_L_bw::PBW
    cr::CT
end

function CollinearPhaseMatch(
    cr::NonlinearCrystal,
    lambda_r1_r2_b::AbstractVector{<:Length},
    temp::Temperature,
    hi_or_lo_r1_r2_b::AbstractVector{<:Symbol},
    theta_pm::Angle,
    phi_pm::Angle,
)
    data_hi_lo_r1 = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[1], temp)
    data_hi_lo_r2 = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[2], temp)
    data_hi_lo_b = refraction_data_hi_lo(theta_pm, phi_pm, cr, lambda_r1_r2_b[3], temp)
    data_hi_lo = [data_hi_lo_r1, data_hi_lo_r2, data_hi_lo_b]

    select_hi_lo = data -> [(hi_or_lo_r1_r2_b[i] == :hi ? data[i][1] : data[i][2]) for i in eachindex(lambda_r1_r2_b)]
    refractive_index_r1_r2_b = select_hi_lo([d[1] for d in data_hi_lo])
    D_dir_r1_r2_b = select_hi_lo([d[2] for d in data_hi_lo])
    E_dir_r1_r2_b = select_hi_lo([d[3] for d in data_hi_lo])
    S_dir_r1_r2_b = select_hi_lo([d[4] for d in data_hi_lo])
    walkoff_angle_r1_r2_b = select_hi_lo([d[5] for d in data_hi_lo])
    group_index_r1_r2_b = select_hi_lo([d[6] for d in data_hi_lo])
    # beta_0_r1_r2_b = select_hi_lo([d[7] for d in data_hi_lo])
    # beta_1_r1_r2_b = select_hi_lo([d[8] for d in data_hi_lo])
    beta_2_r1_r2_b = select_hi_lo([d[9] for d in data_hi_lo])
    beta_3_r1_r2_b = select_hi_lo([d[10] for d in data_hi_lo])

    # Assert that the signs of all field direction vectors (until here chosen randomly) are unified, e.g. all high-index directions point in the same direction 
    k = angles_to_vector(theta_pm, phi_pm)
    first_optical_axis_b = angles_to_vector(optical_axis_angle(cr, lambda_r1_r2_b[3], temp), 0.0u"°")
    unified_dir_signs = calc_unified_dir_signs(E_dir_r1_r2_b, k, first_optical_axis_b)
    D_dir_r1_r2_b = unified_dir_signs .* D_dir_r1_r2_b
    E_dir_r1_r2_b = unified_dir_signs .* E_dir_r1_r2_b

    o_or_e_r1_r2_b = assign_e_o(E_dir_r1_r2_b, cr, use_o_e_for_biaxial=true)

    d_eff = calc_d_eff(cr, E_dir_r1_r2_b...)
    S_0 = (ε_0 * c_0 * lambda_r1_r2_b[1] * lambda_r1_r2_b[2] * prod(refractive_index_r1_r2_b)) / (8 * π^2 * d_eff^2)

    lambda_L_bw = lambda_L_bandwidths(group_index_r1_r2_b)
    temp_L_bw = temperature_L_bandwidth(theta_pm, phi_pm, hi_or_lo_r1_r2_b, cr; lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3], temp)
    theta_L_bw = theta_L_bandwidth(theta_pm, phi_pm, hi_or_lo_r1_r2_b, cr; lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3], temp)
    phi_L_bw = phi_L_bandwidth(theta_pm, phi_pm, hi_or_lo_r1_r2_b, cr; lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3], temp)

    return CollinearPhaseMatch(
        lambda_r1_r2_b,
        temp,
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
        S_0,
        lambda_L_bw,
        temp_L_bw,
        theta_L_bw,
        phi_L_bw,
        cr,
    )
end

function calc_unified_dir_signs(
    dir_vecs_r1_r2_b::AbstractVector{<:AbstractVector{<:Number}},
    k::AbstractVector{<:Number},
    first_optical_axis_b::AbstractVector{<:Number}
)
    # TODO: Check this and the signs of d_eff at opposite hemispheres
    # Use the plane defined by k and the first optical axis of the blue beam to define the blue polarization direction sign; this is supposed to prevent sign jumps, e.g. in d_eff
    dir_ref = cross(first_optical_axis_b, k)
    sign_b = sign(dot(dir_vecs_r1_r2_b[3], dir_ref))

    d_r1_ref = dot(dir_ref, dir_vecs_r1_r2_b[1])
    d_r1_b = dot(dir_vecs_r1_r2_b[3], dir_vecs_r1_r2_b[1])
    d_r2_ref = dot(dir_ref, dir_vecs_r1_r2_b[2])
    d_r2_b = dot(dir_vecs_r1_r2_b[3], dir_vecs_r1_r2_b[2])

    sign_r1 = (abs(d_r1_ref) > abs(d_r1_b)) ? sign(d_r1_ref) : sign(d_r1_b)
    sign_r2 = (abs(d_r2_ref) > abs(d_r2_b)) ? sign(d_r2_ref) : sign(d_r2_b)
    return [sign_r1, sign_r2, sign_b]
end

function assign_e_o(
    E_dir_r1_r2_b::AbstractVector{<:AbstractVector{<:Number}},
    cr::NonlinearCrystal;
    use_o_e_for_biaxial::Bool=true,
    angle_tol=0.2u"°",
)
    o_or_e_r1_r2_b = [nothing, nothing, nothing]
    if cr isa UnidirectionalCrystal
        is_o = [isapprox(acos(clamp(abs(dot([0, 0, 1], E_dir)), 0, 1)), 90.0u"°", atol=ustrip(u"rad", angle_tol)) for E_dir in E_dir_r1_r2_b]
        o_or_e_r1_r2_b = [i ? (:o) : (:e) for i in is_o]
    elseif use_o_e_for_biaxial # Use old notation, where ordinary beams in `BidirectionalCrystal`s are those parallel to one of the principal axes
        is_x_par = [isapprox(acos(clamp(abs(dot([1, 0, 0], E_dir)), 0, 1)), 0.0u"°", atol=ustrip(u"rad", angle_tol)) for E_dir in E_dir_r1_r2_b]
        is_y_par = [isapprox(acos(clamp(abs(dot([0, 1, 0], E_dir)), 0, 1)), 0.0u"°", atol=ustrip(u"rad", angle_tol)) for E_dir in E_dir_r1_r2_b]
        is_z_par = [isapprox(acos(clamp(abs(dot([0, 0, 1], E_dir)), 0, 1)), 0.0u"°", atol=ustrip(u"rad", angle_tol)) for E_dir in E_dir_r1_r2_b]
        o_or_e_r1_r2_b_tmp = [(is_x_par[i] || is_y_par[i] || is_z_par[i]) ? (:o) : (:e) for i in eachindex(E_dir_r1_r2_b)]
        if (:o) in o_or_e_r1_r2_b_tmp
            o_or_e_r1_r2_b = o_or_e_r1_r2_b_tmp # Only set if there is at least on ordinary beam 
        end
    end
    return o_or_e_r1_r2_b
end

function Base.show(io::IO, cpm::CollinearPhaseMatch)
    digits = 3
    println(io, "Crystal: $(cpm.cr.metadata[:description])")
    if !all(isnothing.(cpm.o_or_e_r1_r2_b))
        println(io, "Wavelengths: λ_r1: $(round(u"nm", cpm.lambda_r1_r2_b[1]; digits)) ($(cpm.hi_or_lo_r1_r2_b[1])/$(cpm.o_or_e_r1_r2_b[1])) + λ_r2: $(round(u"nm", cpm.lambda_r1_r2_b[2]; digits)) ($(cpm.hi_or_lo_r1_r2_b[2])/$(cpm.o_or_e_r1_r2_b[2])) = λ_b: $(round(u"nm", cpm.lambda_r1_r2_b[3]; digits)) ($(cpm.hi_or_lo_r1_r2_b[3])/$(cpm.o_or_e_r1_r2_b[3]))")
    else
        println(io, "Wavelengths: λ_r1: $(round(u"nm", cpm.lambda_r1_r2_b[1]; digits)) ($(cpm.hi_or_lo_r1_r2_b[1])) + λ_r2: $(round(u"nm", cpm.lambda_r1_r2_b[2]; digits)) ($(cpm.hi_or_lo_r1_r2_b[2])) = λ_b: $(round(u"nm", cpm.lambda_r1_r2_b[3]; digits)) ($(cpm.hi_or_lo_r1_r2_b[3]))")
    end
    println(io, "Temperature: $(cpm.temp |> u"K") ($(float(cpm.temp |> u"°C")))")
    println(io, "θ: $(round(u"°", cpm.theta_pm |> u"°"; digits)), ϕ: $(round(u"°", cpm.phi_pm |> u"°"; digits))")
    println(io, "Walkoff: $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[1] |> u"mrad"; digits)), $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[2] |> u"mrad"; digits)), $(round(u"mrad", cpm.walkoff_angle_r1_r2_b[3] |> u"mrad"; digits))")
    println(io, "Refractive index: $(round(cpm.refractive_index_r1_r2_b[1]; digits)), $(round(cpm.refractive_index_r1_r2_b[2]; digits)), $(round(cpm.refractive_index_r1_r2_b[3]; digits))")
    println(io, "Group index: $(round(cpm.group_index_r1_r2_b[1]; digits)), $(round(cpm.group_index_r1_r2_b[2]; digits)), $(round(cpm.group_index_r1_r2_b[3]; digits))")
    println(io, "GDD: $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[1]; digits)), $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[2]; digits)), $(round(u"fs^2/mm", cpm.beta_2_r1_r2_b[3]; digits))")
    println(io, "TOD: $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[1]; digits)), $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[2]; digits)), $(round(u"fs^3/mm", cpm.beta_3_r1_r2_b[3]; digits))")
    println(io, "d_eff: $(round(u"pm/V", cpm.d_eff; digits))")
    println(io, "S_0 * L^2: $(round(u"W", cpm.S_0; sigdigits=digits))")
    println(io, "ω_bandwidth(λ_r1=const.) * L: $(round(u"GHz * cm", cpm.lambda_L_bw[1]; digits))")
    println(io, "ω_bandwidth(λ_r2=const.) * L: $(round(u"GHz * cm", cpm.lambda_L_bw[2]; digits))")
    println(io, "ω_bandwidth(λ_b=const.) * L: $(round(u"GHz * cm", cpm.lambda_L_bw[3]; digits))")
    println(io, "T_bandwidth * L: $(round(u"K * cm", cpm.temp_L_bw; digits))")
    println(io, "θ_bandwidth * L: $(round(u"mrad * cm", cpm.theta_L_bw; digits))")
    println(io, "ϕ_bandwidth * L: $(round(u"mrad * cm", cpm.phi_L_bw; digits))")
    println(io, "E_dir: $(round.(cpm.E_dir_r1_r2_b[1]; digits)), $(round.(cpm.E_dir_r1_r2_b[2]; digits)), $(round.(cpm.E_dir_r1_r2_b[3]; digits))")
    # println(io, "D_dir: $(round.(cpm.D_dir_r1_r2_b[1]; digits)), $(round.(cpm.D_dir_r1_r2_b[2]; digits)), $(round.(cpm.D_dir_r1_r2_b[3]; digits))")
    # println(io, "E_dir_angles: $(round.(u"°", vector_to_angles(cpm.E_dir_r1_r2_b[1]); digits)), $(round.(u"°", vector_to_angles(cpm.E_dir_r1_r2_b[2]); digits)), $(round.(u"°", vector_to_angles(cpm.E_dir_r1_r2_b[3]); digits))")
    # println(io, "D_dir_angles: $(round.(u"°", vector_to_angles(cpm.D_dir_r1_r2_b[1]); digits)), $(round.(u"°", vector_to_angles(cpm.D_dir_r1_r2_b[2]); digits)), $(round.(u"°", vector_to_angles(cpm.D_dir_r1_r2_b[3]); digits))")
end

function pm_wavelengths(;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
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

    @assert 1 / lambda_r1 + 1 / lambda_r2 ≈ 1 / lambda_b "Phasematching is only fulfilled for 1 / λ_r1 + 1 / λ_r2 = 1 / λ_b."

    return lambda_r1, lambda_r2, lambda_b
end

function delta_k(
    θ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::UnidirectionalCrystal;
    kwargs...
)
    # Use symmetry in unidirectional crystals
    return delta_k(θ_pm, 0.0u"°", hi_or_lo_r1_r2_b, cr; kwargs...)
end

function delta_k(
    θ_pm::Angle,
    ϕ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_r1_r2_b])

    n_r1 = hi_or_lo_r1_r2_b[1] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, temp; n_hi_lo_only=true)[2]
    n_r2 = hi_or_lo_r1_r2_b[2] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, temp; n_hi_lo_only=true)[2]
    n_b = hi_or_lo_r1_r2_b[3] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, temp; n_hi_lo_only=true)[2]

    return 2π * (n_r1 / lambda_r1 + n_r2 / lambda_r2 - n_b / lambda_b)
end

function delta_k_with_shifting(
    θ_pm::Angle,
    ϕ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref,
    delta_theta=0.0u"°",
    delta_phi=0.0u"°",
    delta_temp::Temperature=0.0u"K",
    delta_omega_r1::Frequency=0.0u"GHz",
    delta_omega_r2::Frequency=0.0u"GHz",
    delta_omega_b::Frequency=0.0u"GHz",
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_r1_r2_b])

    @assert delta_omega_r1 + delta_omega_r2 == delta_omega_b "Phasematching is only fulfilled for Δω_r1 + Δω_r2 = Δω_b"
    lambda_r1_shifted = shift_lambda_with_freq(lambda_r1, delta_omega_r1)
    lambda_r2_shifted = shift_lambda_with_freq(lambda_r2, delta_omega_r2)
    lambda_b_shifted = shift_lambda_with_freq(lambda_b, delta_omega_b)

    delta_k(
        (θ_pm + delta_theta) |> u"rad",
        (ϕ_pm + delta_phi) |> u"rad",
        hi_or_lo_r1_r2_b,
        cr;
        lambda_r1=lambda_r1_shifted,
        lambda_r2=lambda_r2_shifted,
        lambda_b=lambda_b_shifted,
        temp=temp + delta_temp,
    )
end

# function delta_k(
#     cpm::CollinearPhaseMatch;
#     delta_theta=0.0u"°",
#     delta_phi=0.0u"°",
#     delta_temp::Temperature=0.0u"K",
#     delta_omega_r1::Frequency=0.0u"GHz",
#     delta_omega_r2::Frequency=0.0u"GHz",
#     delta_omega_b::Frequency=0.0u"GHz",
# )
#     @assert delta_omega_r1 + delta_omega_r2 == delta_omega_b "Phasematching is only fulfilled for Δω_r1 + Δω_r2 = Δω_b"
#     lambda_r1_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[1], delta_omega_r1)
#     lambda_r2_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[2], delta_omega_r2)
#     lambda_b_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[3], delta_omega_b)

#     return delta_k(
#         (cpm.theta_pm + delta_theta) |> u"rad",
#         (cpm.phi_pm + delta_phi) |> u"rad",
#         cpm.hi_or_lo_r1_r2_b,
#         cpm.cr;
#         lambda_r1=lambda_r1_shifted,
#         lambda_r2=lambda_r2_shifted,
#         lambda_b=lambda_b_shifted,
#         temp=cpm.temp + delta_temp,
#     )
# end

function lambda_L_bandwidths(group_index_r1_r2_b::AbstractVector{<:Real})
    # ∂Δk/∂f_r2 = ∂k_r2/∂f_r2 - ∂k_b/∂f_b for f_r1=const.
    ddelta_k_df_r2_neg_b = 2π / c_0 * (group_index_r1_r2_b[2] - group_index_r1_r2_b[3])

    # ∂Δk/∂f_r1 = ∂k_r1/∂f_r1 - ∂k_b/∂f_b for f_r2=const.
    ddelta_k_df_r1_neg_b = 2π / c_0 * (group_index_r1_r2_b[1] - group_index_r1_r2_b[3])

    # ∂Δk/∂f_r1 = ∂k_r1/∂f_r1 - ∂k_r2/∂f_r2 for f_b=const.
    ddelta_k_df_r1_neg_r2 = 2π / c_0 * (group_index_r1_r2_b[1] - group_index_r1_r2_b[2])

    # Δf_k_bw  * L = 2π / |∂k_i/∂f_i - ∂k_j/∂f_j| = 2π / |1 / v_g_i - 1 / v_g_j|
    return (
        2π / abs(ddelta_k_df_r2_neg_b),
        2π / abs(ddelta_k_df_r1_neg_b),
        2π / abs(ddelta_k_df_r1_neg_r2)
    )
end


function temperature_L_bandwidth(
    θ_pm::Angle,
    ϕ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref
)
    fun = ΔT -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            θ_pm,
            ϕ_pm,
            hi_or_lo_r1_r2_b,
            cr;
            lambda_r1,
            lambda_r2,
            lambda_b,
            temp,
            delta_temp=ΔT * u"K",
        )
    )
    ddeltak_dT = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * K^-1"

    # ΔT_bw * L = 2π / |∂Δk/∂T|
    return 2π / abs(ddeltak_dT)
end

function theta_L_bandwidth(
    θ_pm::Angle,
    ϕ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref
)
    fun = Δθ -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            θ_pm,
            ϕ_pm,
            hi_or_lo_r1_r2_b,
            cr;
            lambda_r1,
            lambda_r2,
            lambda_b,
            temp,
            delta_theta=Δθ * u"rad",
        )
    )
    ddeltak_dtheta = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * rad^-1"

    # Δθ_bw * L = 2π / |∂Δk/∂θ|
    return 2π / abs(ddeltak_dtheta)
end

function phi_L_bandwidth(
    θ_pm::Angle,
    ϕ_pm::Angle,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref
)
    fun = Δϕ -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            θ_pm,
            ϕ_pm,
            hi_or_lo_r1_r2_b,
            cr;
            lambda_r1,
            lambda_r2,
            lambda_b,
            temp,
            delta_phi=Δϕ * u"rad",
        )
    )
    ddeltak_dphi = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * rad^-1"

    # Δϕ_bw * L = 2π / |∂Δk/∂ϕ|
    return 2π / abs(ddeltak_dphi)
end


# function find_pms(
#     pol_r1_r2_b::AbstractVector{Symbol},
#     cr::NonlinearCrystal;
#     lambda_r1::Union{Nothing,Length}=nothing,
#     lambda_r2::Union{Nothing,Length}=nothing,
#     lambda_b::Union{Nothing,Length}=nothing,
#     temp::Temperature=default_temp(cr),
#     theta_fixed=nothing,
#     phi_fixed=nothing,
#     ngrid=500,
#     tol=1e-14u"nm^-1",
# )
#     lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
#     lambda_r1_r2_b = [lambda_r1, lambda_r2, lambda_b]

#     hi_or_lo_r1_r2_b = polarization_r1_r2_b_to_hi_lo(pol_r1_r2_b, cr, lambda_r1_r2_b; temp)

#     pms = CollinearPhaseMatch[]

#     if !isnothing(theta_fixed) && isnothing(phi_fixed)
#         # Global search over ϕ
#         all_ϕ = range(0, 2π, length=ngrid) * u"rad"
#         all_delta_k = [delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for ϕ in all_ϕ]

#         # Find approximate zero-crossings
#         for i in 1:(length(all_ϕ)-1)
#             if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0  # zero crossing detected
#                 # Refine solution using local optimization
#                 ϕ_sol = find_zero(ϕ -> delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
#                     (all_ϕ[i], all_ϕ[i+1]), Bisection(), atol=tol)
#                 push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, theta_fixed, ϕ_sol))
#             end
#         end

#     elseif isnothing(theta_fixed) && !isnothing(phi_fixed)
#         all_θ = range(0, π, length=ngrid) * u"rad"
#         all_delta_k = [delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in all_θ]

#         # Find approximate zero-crossings
#         for i in 1:(length(all_θ)-1)
#             if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0 # zero crossing detected
#                 # Refine solution using local optimization
#                 θ_sol = find_zero(θ -> delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
#                     (all_θ[i], all_θ[i+1]), Bisection(), atol=tol)
#                 push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, θ_sol, phi_fixed))
#             end
#         end
#     else
#         error("You must provide either theta_fixed or phi_fixed.")
#     end
#     return pms
# end

# function find_pms(
#     cr::NonlinearCrystal;
#     lambda_r1::Union{Nothing,Length}=nothing,
#     lambda_r2::Union{Nothing,Length}=nothing,
#     lambda_b::Union{Nothing,Length}=nothing,
#     temp::Temperature=default_temp(cr),
#     pol_r1_r2_b_vec::Union{Nothing,AbstractVector{<:AbstractVector{Symbol}}}=nothing,
#     theta_fixed=nothing,
#     phi_fixed=nothing,
#     ngrid=500,
#     tol=1e-14u"nm^-1",
# )
#     if isnothing(pol_r1_r2_b_vec)
#         pol_r1_r2_b_vec = bool_permutations(:hi, :lo, 3)
#     end

#     all_res = map(pol_r1_r2_b_vec) do hi_or_lo_r1_r2_b
#         find_pms(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp, theta_fixed, phi_fixed, ngrid, tol)
#     end

#     return reduce(vcat, all_res)
# end

function find_nearest_pm_along_theta_phi(
    theta_target::Angle,
    phi_target::Angle,
    pol_r1_r2_b_vec::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    all_pm_θ_fixed = find_all_pms_along_dimension(pol_r1_r2_b_vec, cr; lambda_r1_fixed=lambda_r1, lambda_r2_fixed=lambda_r2, lambda_b_fixed=lambda_b, temp_min=temp, temp_max=temp, ngrid, tol, theta_fixed=theta_target)
    all_pm_ϕ_fixed = find_all_pms_along_dimension(pol_r1_r2_b_vec, cr; lambda_r1_fixed=lambda_r1, lambda_r2_fixed=lambda_r2, lambda_b_fixed=lambda_b, temp_min=temp, temp_max=temp, ngrid, tol, phi_fixed=phi_target)

    all_pm_candidates = [[pm for pm in all_pm_θ_fixed]; [pm for pm in all_pm_ϕ_fixed]]
    isempty(all_pm_candidates) && return nothing

    all_θ_pm_candidates = [pm.theta_pm for pm in all_pm_candidates]
    all_ϕ_pm_candidates = [pm.phi_pm for pm in all_pm_candidates]

    target_vec = angles_to_vector(theta_target, phi_target)
    pm_candidate_vecs = angles_to_vector.(all_θ_pm_candidates, all_ϕ_pm_candidates)
    i_nearest = findmax([dot(pm_vec, target_vec) for pm_vec in pm_candidate_vecs])[2]

    return all_pm_candidates[i_nearest]
end

function find_nearest_pm_along_lambda_r_b(
    pol_r1_r2_b_vec::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    pricipal_axes::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    @assert count(isnothing.([lambda_r1, lambda_r2])) == 1 "Either lambda_r1 or lambda_r2 must be given and the other one must be nothing."

    # Search along lambda_b
    all_pm_lambda_b = find_all_pms_along_dimension(pol_r1_r2_b_vec, cr;
        lambda_b_fixed=lambda_b,
        temp_min=temp,
        temp_max=temp,
        pricipal_axes=pricipal_axes,
        ngrid=ngrid,
        tol=tol)

    # Search along lambda_r1 or lambda_r2
    all_pm_lambda_r = find_all_pms_along_dimension(pol_r1_r2_b_vec, cr;
        lambda_r1_fixed=lambda_r1,
        lambda_r2_fixed=lambda_r2,
        temp_min=temp,
        temp_max=temp,
        pricipal_axes=pricipal_axes,
        ngrid=ngrid,
        tol=tol)

    # Combine candidates
    all_pm_candidates = vcat(all_pm_lambda_b, all_pm_lambda_r)
    isempty(all_pm_candidates) && return nothing

    # Compute distance to target wavelength(s)
    distances = map(all_pm_candidates) do pm
        λs = pm.lambda_r1_r2_b
        λ_r = isnothing(lambda_r1) ? λs[2] : λs[1]  # whichever was scanned
        λ_b = λs[3]
        Δλ_r = abs(λ_r - (lambda_r1 === nothing ? lambda_r2 : lambda_r1))
        Δλ_b = abs(λ_b - lambda_b)
        Δλ_r^2 + Δλ_b^2
    end

    i_nearest = findmin(distances)[2]
    return all_pm_candidates[i_nearest]
end

function find_all_pms_along_dimension(
    pol_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1_fixed::Union{Nothing,Length}=nothing,
    lambda_r2_fixed::Union{Nothing,Length}=nothing,
    lambda_b_fixed::Union{Nothing,Length}=nothing,
    temp_min::Temperature=default_temp(cr),
    temp_max::Temperature=default_temp(cr),
    pricipal_axes::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    theta_fixed=nothing,
    phi_fixed=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    λmin, λmax = valid_lambda_range(cr)
    temp_fixed = (temp_min == temp_max) ? temp_min : nothing

    if !isnothing(theta_fixed) || !isnothing(phi_fixed)
        @assert isnothing(pricipal_axes) "pricipal_axes must be `nothing` if searching along fixed θ or ϕ directions."
        @assert !isnothing(temp_fixed) "Fix temperature by selecting temp_min = temp_max if searching along fixed θ or ϕ directions."
        return _pms_vs_angle(pol_r1_r2_b, cr; lambda_r1=lambda_r1_fixed, lambda_r2=lambda_r2_fixed, lambda_b=lambda_b_fixed, temp=temp_fixed, theta_fixed, phi_fixed, ngrid, tol)
    end
    
    red_fixed_count = count(.!isnothing.([lambda_r1_fixed, lambda_r2_fixed]))
    @assert red_fixed_count < 2 "Only either lambda_r1_fixed or lambda_r2_fixed should be provided but both are given."
    @assert (count(.!isnothing.([temp_fixed, lambda_b_fixed])) + red_fixed_count) == 2 "You must provide exactly two of those: (temp_min == temp_max), lambda_b_fixed, and (lambda_r1_fixed or lambda_r2_fixed)."

    hi_or_lo_r1_r2_b = pol_r1_r2_b # TODO: Enable searching along temperature and wavelengths with :o/:e notation

    if isnothing(pricipal_axes) && isnothing(theta_fixed) && isnothing(phi_fixed)
        pricipal_axes = [:X, :Y, :Z]
    end
    θ_ϕ_pricipal_axes = axes_to_θ_ϕ(pricipal_axes)

    if isnothing(temp_fixed)
        return _pms_vs_temperature(θ_ϕ_pricipal_axes, hi_or_lo_r1_r2_b, cr, lambda_r1_fixed, lambda_r2_fixed, lambda_b_fixed, temp_min, temp_max, ngrid, tol)
    elseif isnothing(lambda_b_fixed)
        return _pms_vs_lambda_b(θ_ϕ_pricipal_axes, hi_or_lo_r1_r2_b, cr, lambda_r1_fixed, lambda_r2_fixed, temp_fixed, λmin, λmax, ngrid, tol)
    elseif isnothing(lambda_r1_fixed) && isnothing(lambda_r2_fixed)
        return _pms_vs_lambda_r1(θ_ϕ_pricipal_axes, hi_or_lo_r1_r2_b, cr, lambda_b_fixed, temp_fixed, λmin, λmax, ngrid, tol)
    else
        error("This should never be reached.")
    end
end

function _pms_vs_angle(pol_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp, theta_fixed, phi_fixed, ngrid, tol)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    lambda_r1_r2_b = [lambda_r1, lambda_r2, lambda_b]

    # Transform :o/e into :hi/:lo if possible
    hi_or_lo_r1_r2_b = polarization_r1_r2_b_to_hi_lo(pol_r1_r2_b, cr, lambda_r1_r2_b; temp)

    pms = CollinearPhaseMatch[]

    if !isnothing(theta_fixed) && isnothing(phi_fixed)
        all_ϕ = range(0, 2π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for ϕ in all_ϕ]

        for i in 1:(length(all_ϕ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0
                ϕ_sol = find_zero(ϕ -> delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_ϕ[i], all_ϕ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, theta_fixed, ϕ_sol))
            end
        end

    elseif isnothing(theta_fixed) && !isnothing(phi_fixed)
        all_θ = range(0, π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in all_θ]

        for i in 1:(length(all_θ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0
                θ_sol = find_zero(θ -> delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_θ[i], all_θ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, θ_sol, phi_fixed))
            end
        end
    else
        error("You must provide either theta_fixed or phi_fixed (but not both).")
    end

    return pms
end

function _pms_vs_temperature(θ_ϕs, hi_or_lo_r1_r2_b, cr, λr1_fix, λr2_fix, λb_fix, temp_min, temp_max, ngrid, tol)
    pms = CollinearPhaseMatch[]
    λr1, λr2, λb = pm_wavelengths(; lambda_r1=λr1_fix, lambda_r2=λr2_fix, lambda_b=λb_fix)

    for θ_ϕ in θ_ϕs
        temps = range(temp_min, temp_max, length=ngrid)
        Δk = [delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; temp=T, lambda_r1=λr1, lambda_r2=λr2, lambda_b=λb) for T in temps]

        for i in 1:length(temps)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                fun = T -> delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; temp=T, lambda_r1=λr1, lambda_r2=λr2, lambda_b=λb)
                T_sol = find_zero(fun, (temps[i], temps[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, [λr1, λr2, λb], T_sol, hi_or_lo_r1_r2_b, θ_ϕ...))
            end
        end
    end

    return pms
end

function _pms_vs_lambda_b(θ_ϕs, hi_or_lo_r1_r2_b, cr, λr1_fix, λr2_fix, temp, λmin, λmax, ngrid, tol)
    pms = CollinearPhaseMatch[]

    for θ_ϕ in θ_ϕs
        if !isnothing(λr1_fix)
            λb_range = range(
                max(λmin, 1 / (1 / λr1_fix + 1 / λmin)),
                min(λmax, 1 / (1 / λr1_fix + 1 / λmax)),
                length=ngrid
            )
            Δk = [delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=λr1_fix, lambda_r2=(1 / (1 / λb - 1 / λr1_fix)), lambda_b=λb, temp=temp) for λb in λb_range]
        else
            λb_range = range(
                max(λmin, 1 / (1 / λr2_fix + 1 / λmin)),
                min(λmax, 1 / (1 / λr2_fix + 1 / λmax)),
                length=ngrid
            )
            Δk = [delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=(1 / (1 / λb - 1 / λr2_fix)), lambda_r2=λr2_fix, lambda_b=λb, temp=temp) for λb in λb_range]
        end

        for i in 1:length(λb_range)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                if !isnothing(λr1_fix)
                    fun = λb -> delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=λr1_fix, lambda_r2=(1 / (1 / λb - 1 / λr1_fix)), lambda_b=λb, temp=temp)
                    λb_sol = find_zero(fun, (λb_range[i], λb_range[i+1]), Bisection(), atol=tol)
                    λr2_sol = 1 / (1 / λb_sol - 1 / λr1_fix)
                    push!(pms, CollinearPhaseMatch(cr, [λr1_fix, λr2_sol, λb_sol], temp, hi_or_lo_r1_r2_b, θ_ϕ...))
                else
                    fun = λb -> delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=(1 / (1 / λb - 1 / λr2_fix)), lambda_r2=λr2_fix, lambda_b=λb, temp=temp)
                    λb_sol = find_zero(fun, (λb_range[i], λb_range[i+1]), Bisection(), atol=tol)
                    λr1_sol = 1 / (1 / λb_sol - 1 / λr2_fix)
                    push!(pms, CollinearPhaseMatch(cr, [λr1_sol, λr2_fix, λb_sol], temp, hi_or_lo_r1_r2_b, θ_ϕ...))
                end
            end
        end
    end

    return pms
end

function _pms_vs_lambda_r1(θ_ϕs, hi_or_lo_r1_r2_b, cr, λb, temp, λmin, λmax, ngrid, tol)
    pms = CollinearPhaseMatch[]

    for θ_ϕ in θ_ϕs
        λr1_range = range(
            max(λmin, 1 / (1 / λb - 1 / λmax)),
            min(2 * λb, 1 / (1 / λb - 1 / (2 * λb))),
            length=ngrid
        )
        Δk = [delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=λr1, lambda_r2=(1 / (1 / λb - 1 / λr1)), lambda_b=λb, temp=temp) for λr1 in λr1_range]

        for i in 1:length(λr1_range)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                fun = λr1 -> delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1=λr1, lambda_r2=(1 / (1 / λb - 1 / λr1)), lambda_b=λb, temp=temp)
                λr1_sol = find_zero(fun, (λr1_range[i], λr1_range[i+1]), Bisection(), atol=tol)
                λr2_sol = 1 / (1 / λb - 1 / λr1_sol)
                push!(pms, CollinearPhaseMatch(cr, [λr1_sol, λr2_sol, λb], temp, hi_or_lo_r1_r2_b, θ_ϕ...))
            end
        end
    end

    return pms
end