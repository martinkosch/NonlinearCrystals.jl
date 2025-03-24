export PhaseMatch, CollinearPhaseMatch, pm_wavelengths, find_noncritical_pms, find_pms, find_nearest_pm, delta_k, plot_delta_k_map, plot_noncritical_pms

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
    theta_pm,
    phi_pm,
)
    theta_pm = theta_pm |> u"rad"
    phi_pm = phi_pm |> u"rad"

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
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=cr.n_x_principal.temp_ref
)
    θ_pm = θ_pm |> u"rad"
    ϕ_pm = ϕ_pm |> u"rad"
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_r1_r2_b])

    n_r1 = hi_or_lo_r1_r2_b[1] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r1, temp; n_hi_lo_only=true)[2]
    n_r2 = hi_or_lo_r1_r2_b[2] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_r2, temp; n_hi_lo_only=true)[2]
    n_b = hi_or_lo_r1_r2_b[3] == :hi ? refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, temp; n_hi_lo_only=true)[1] : refraction_data_hi_lo(θ_pm, ϕ_pm, cr, lambda_b, temp; n_hi_lo_only=true)[2]

    return 2π * (n_r1 / lambda_r1 + n_r2 / lambda_r2 - n_b / lambda_b)
end

function delta_k_with_shifting(
    θ_pm,
    ϕ_pm,
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
    θ_pm_shifted = (θ_pm + delta_theta) |> u"rad"
    ϕ_pm_shifted = (ϕ_pm + delta_phi) |> u"rad"

    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_r1_r2_b])

    @assert delta_omega_r1 + delta_omega_r2 == delta_omega_b "Phasematching is only fulfilled for Δω_r1 + Δω_r2 = Δω_b"
    lambda_r1_shifted = shift_lambda_with_freq(lambda_r1, delta_omega_r1)
    lambda_r2_shifted = shift_lambda_with_freq(lambda_r2, delta_omega_r2)
    lambda_b_shifted = shift_lambda_with_freq(lambda_b, delta_omega_b)

    delta_k(
        θ_pm_shifted,
        ϕ_pm_shifted,
        hi_or_lo_r1_r2_b,
        cr;
        lambda_r1=lambda_r1_shifted,
        lambda_r2=lambda_r2_shifted,
        lambda_b=lambda_b_shifted,
        temp=temp + delta_temp,
    )
end

function delta_k(
    cpm::CollinearPhaseMatch;
    delta_theta=0.0u"°",
    delta_phi=0.0u"°",
    delta_temp::Temperature=0.0u"K",
    delta_omega_r1::Frequency=0.0u"GHz",
    delta_omega_r2::Frequency=0.0u"GHz",
    delta_omega_b::Frequency=0.0u"GHz",
)
    @assert delta_omega_r1 + delta_omega_r2 == delta_omega_b "Phasematching is only fulfilled for Δω_r1 + Δω_r2 = Δω_b"
    lambda_r1_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[1], delta_omega_r1)
    lambda_r2_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[2], delta_omega_r2)
    lambda_b_shifted = shift_lambda_with_freq(cpm.lambda_r1_r2_b[3], delta_omega_b)

    return delta_k(
        cpm.theta_pm + delta_theta,
        cpm.phi_pm + delta_phi,
        cpm.hi_or_lo_r1_r2_b,
        cpm.cr;
        lambda_r1=lambda_r1_shifted,
        lambda_r2=lambda_r2_shifted,
        lambda_b=lambda_b_shifted,
        temp=cpm.temp + delta_temp,
    )
end

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
    θ_pm,
    ϕ_pm,
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
    θ_pm,
    ϕ_pm,
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
    θ_pm,
    ϕ_pm,
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


function find_pms(
    pol_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    theta_fixed=nothing,
    phi_fixed=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    lambda_r1_r2_b = [lambda_r1, lambda_r2, lambda_b]

    hi_or_lo_r1_r2_b = polarization_r1_r2_b_to_hi_lo(pol_r1_r2_b, cr, lambda_r1_r2_b; temp)

    pms = CollinearPhaseMatch[]

    if theta_fixed !== nothing && isnothing(phi_fixed)
        # Global search over ϕ
        all_ϕ = range(0, 2π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for ϕ in all_ϕ]

        # Find approximate zero-crossings
        for i in 1:(length(all_ϕ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0  # zero crossing detected
                # Refine solution using local optimization
                ϕ_sol = find_zero(ϕ -> delta_k(theta_fixed, ϕ, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_ϕ[i], all_ϕ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, theta_fixed, ϕ_sol))
            end
        end

    elseif isnothing(theta_fixed) && !isnothing(phi_fixed)
        all_θ = range(0, π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in all_θ]

        # Find approximate zero-crossings
        for i in 1:(length(all_θ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0 # zero crossing detected
                # Refine solution using local optimization
                θ_sol = find_zero(θ -> delta_k(θ, phi_fixed, hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_θ[i], all_θ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, temp, hi_or_lo_r1_r2_b, θ_sol, phi_fixed))
            end
        end
    else
        error("You must provide either theta_fixed or phi_fixed.")
    end
    return pms
end

function find_pms(
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    pol_r1_r2_b_vec::Union{Nothing,AbstractVector{<:AbstractVector{Symbol}}}=nothing,
    theta_fixed=nothing,
    phi_fixed=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    if isnothing(pol_r1_r2_b_vec)
        pol_r1_r2_b_vec = bool_permutations(:hi, :lo, 3)
    end

    all_res = map(pol_r1_r2_b_vec) do hi_or_lo_r1_r2_b
        find_pms(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp, theta_fixed, phi_fixed, ngrid, tol)
    end

    return reduce(vcat, all_res)
end

function find_nearest_pm(
    θ_target,
    ϕ_target,
    pol_r1_r2_b_vec::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    θ_target = θ_target |> u"rad"
    ϕ_target = ϕ_target |> u"rad"

    all_pm_θ_fixed = find_pms(pol_r1_r2_b_vec, cr; lambda_r1, lambda_r2, lambda_b, temp, ngrid, tol, theta_fixed=θ_target)
    all_pm_ϕ_fixed = find_pms(pol_r1_r2_b_vec, cr; lambda_r1, lambda_r2, lambda_b, temp, ngrid, tol, phi_fixed=ϕ_target)

    all_pm_candidates = [[pm for pm in all_pm_θ_fixed]; [pm for pm in all_pm_ϕ_fixed]]
    isempty(all_pm_candidates) && return nothing

    all_θ_pm_candidates = [pm.theta_pm for pm in all_pm_candidates]
    all_ϕ_pm_candidates = [pm.phi_pm for pm in all_pm_candidates]

    target_vec = angles_to_vector(θ_target, ϕ_target)
    pm_candidate_vecs = angles_to_vector.(all_θ_pm_candidates, all_ϕ_pm_candidates)
    i_nearest = findmax([dot(pm_vec, target_vec) for pm_vec in pm_candidate_vecs])[2]

    return all_pm_candidates[i_nearest]
end



function find_noncritical_pms(
    pol_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp_min::Temperature=293.15u"K",
    temp_max::Temperature=600u"K",
    axes::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    lambda_r1_r2_b = [lambda_r1, lambda_r2, lambda_b]

    hi_or_lo_r1_r2_b = polarization_r1_r2_b_to_hi_lo(pol_r1_r2_b, cr, lambda_r1_r2_b; temp)

    pms = CollinearPhaseMatch[]

    if isnothing(axes)
        axes = [:X, :Y, :Z]
    end

    # Only search along positive direction of the specified axes
    θ_ϕ_dirs = axes_to_θ_ϕ(axes)
    for θ_ϕ in θ_ϕ_dirs
        # Global search over all temperatures
        all_temp = range(temp_min, temp_max, length=ngrid)
        all_delta_k = [delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b) for temp in all_temp]

        # Find approximate zero-crossings
        for i in 1:(length(all_temp)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0  # zero crossing detected
                # Refine solution using local optimization
                T_sol = find_zero(temp -> delta_k(θ_ϕ..., hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_temp[i], all_temp[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_r1_r2_b, T_sol, hi_or_lo_r1_r2_b, θ_ϕ...))
            end
        end
    end

    return pms
end

function find_noncritical_pms(
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp_min::Temperature=293.15u"K",
    temp_max::Temperature=600u"K",
    axes::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    pol_r1_r2_b_vec::Union{Nothing,AbstractVector{<:AbstractVector{Symbol}}}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    if isnothing(pol_r1_r2_b_vec)
        pol_r1_r2_b_vec = bool_permutations(:hi, :lo, 3)
    end

    all_res = map(pol_r1_r2_b_vec) do hi_or_lo_r1_r2_b
        find_noncritical_pms(hi_or_lo_r1_r2_b, cr; lambda_r1, lambda_r2, lambda_b, temp_min, temp_max, axes, ngrid, tol)
    end

    return reduce(vcat, all_res)
end

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
    pm = find_nearest_pm(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b)
    return isnothing(pm) ? "" : plot_pm_label(pm)
end

function calc_raw_noncritical_pm_lines(
    pricipal_axis::Symbol,
    hi_or_lo_r1_r2_b::AbstractVector{Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
)
    @assert pricipal_axis in [:X, :Y, :Z]

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
            all_delta_k[i_b, i_r12] = ustrip(u"m^-1", delta_k(axes_to_θ_ϕ(pricipal_axis)[1]..., hi_or_lo_r1_r2_b, cr; temp, lambda_r1, lambda_r2, lambda_b))
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
        cb_interp = cb[switch_idx] * (1 - frac) + cb[switch_idx + 1] * frac
        cr_interp = cr[switch_idx] * (1 - frac) + cr[switch_idx + 1] * frac

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

function refined_noncritical_pm_lines(
    pricipal_axis::Symbol,
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
        pricipal_axis, 
        hi_or_lo_r1_r2_b, 
        cr; 
        lambda_b_min, 
        lambda_b_max, 
        lambda_r12_min, 
        lambda_r12_max, 
        temp, 
        ngrid
    )
    
    # If both red waves are polarized equally, the plot lines show a helpful phase match symmetry around the line cont_r_raw = 2 * cont_b_raw
    is_pol_r1_r2_equal = hi_or_lo_r1_r2_b[1] == hi_or_lo_r1_r2_b[1] 
    
    for cl in cpl.lines
        cont_b_raw = [v[1] for v in cl.vertices]
        cont_r_raw = [v[2] for v in cl.vertices]

        for (cb, cr) in zip(split_on_nan(cont_b_raw, cont_r_raw)...)
            shear_cr = cr .- 2 * cb # All SHG points are now on the y = 0 axis
            segments_cb, segments_cr, cb_intersections, cr_intersections, segment_signs = split_and_interpolate(cb, cr, sign_switch_fractions(shear_cr)...)
            
            if is_pol_r1_r2_equal
                # Throw away lower part, as it is symmetric anyway
                
            else
                # Flip lower parts to upper half using phase match symmetry
                
            end

            # @show cb_intersections
            # for i in eachindex(segments_cb)
            #     lines!(segments_cb[i], segments_cr[i] .- 2 * segments_cb[i])
            # end
            # !isempty(cb_intersections) && scatter!(cb_intersections, cr_intersections .- 2 * cb_intersections)
            
            # if !isempty(sign_switches)
            #     scatter!(cb[sign_switches], cr[sign_switches] .- 2 * cb[sign_switches])
            # end
            # lines!(cb, 1 ./ (1 ./ cb - 1 ./ cr))
        end
    end
end

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
    pm = find_nearest_pm(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
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
    pm = find_nearest_pm(θ_pm, ϕ_pm, hi_or_lo_r1_r2_b, cr; temp, lambda_r1=lambda_r1_r2_b[1], lambda_r2=lambda_r1_r2_b[2], lambda_b=lambda_r1_r2_b[3])
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