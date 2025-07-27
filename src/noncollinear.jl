export delta_k_noncollinear, NoncollinearPhaseMatch, to_yaml

struct PMNoncollinearData{CT<:NonlinearCrystal}
    hi_or_lo_rrb::NTuple{3,Symbol}
    pm_type::NTuple{2,Union{Nothing,PMType}}
    theta_pm_rrb::NTuple{3,typeof(1.0u"rad")}
    phi_pm_rrb::NTuple{3,typeof(1.0u"rad")}
    cr::CT
    lambda_rrb::NTuple{3,typeof(1.0u"m")}
    temp::typeof(1.0u"K")
end

function PMNoncollinearData(
    refr_data::PMRefractionData,
    hi_or_lo_rrb::NTuple{3,Symbol},
    theta_pm_rrb::NTuple{3,Angle},
    phi_pm_rrb::NTuple{3,Angle},
    cr::NonlinearCrystal,
    lambda_rrb::NTuple{3,Length},
    temp::Temperature;
    angle_tol::Angle=0.1u"°",
)
    if isa(cr, UnidirectionalCrystal)
        # o and e are valid classifiers for UnidirectionalCrystal even if propagating outside of the pricipal planes
        pm_type = (PMType(:UD, refr_data), nothing)
    else
        principal_planes_rrb = [find_principal_planes(theta_pm_rrb[i], phi_pm_rrb[i]; angle_tol) for i in eachindex(lambda_rrb)]
        @show principal_planes_rrb

        # Only return principal plane if it is identical for all wavelengths
        common(x) = all(x .== x[begin]) ? x[begin] : nothing
        principal_planes = (common([p[i] for p in principal_planes_rrb]) for i in eachindex(principal_planes_rrb[1]))

        pm_type = Tuple((isnothing(p) ? nothing : PMType(p, refr_data) for p in principal_planes))
    end

    return PMNoncollinearData{typeof(cr)}(
        hi_or_lo_rrb,
        pm_type,
        theta_pm_rrb,
        phi_pm_rrb,
        cr,
        lambda_rrb,
        temp,
    )
end

function PMEfficiencyData(pm_data::PMNoncollinearData, refr_data::PMRefractionData)
    d_eff = calc_d_eff(pm_data.cr, refr_data.E_dir_rrb...; lambda_rrb=pm_data.lambda_rrb, temp=pm_data.temp, use_miller_scaling=true)
    d_eff_no_miller = calc_d_eff(pm_data.cr, refr_data.E_dir_rrb...; lambda_rrb=pm_data.lambda_rrb, temp=pm_data.temp, use_miller_scaling=false)
    S_0_Lsquared = (ε_0 * c_0 * pm_data.lambda_rrb[1] * pm_data.lambda_rrb[2] * prod(refr_data.n_rrb)) / (8 * π^2 * d_eff^2)
    S_0_Lsquared_no_miller = (ε_0 * c_0 * pm_data.lambda_rrb[1] * pm_data.lambda_rrb[2] * prod(refr_data.n_rrb)) / (8 * π^2 * d_eff_no_miller^2)
    return PMEfficiencyData(
        d_eff,
        d_eff_no_miller,
        S_0_Lsquared,
        S_0_Lsquared_no_miller,
    )
end

function PMBandwidthData(pm_data::PMNoncollinearData, refr_data::PMRefractionData)
    omega_L_bw = lambda_L_bandwidths(refr_data) # TODO: This is possibly slightly wrong if the formula for collinear phase-matches is used
    temp_L_bw = temperature_L_bandwidth(pm_data)
    theta_L_bw = theta_L_bandwidth(pm_data)
    phi_L_bw = phi_L_bandwidth(pm_data)
    return PMBandwidthData(
        omega_L_bw,
        temp_L_bw,
        theta_L_bw,
        phi_L_bw,
    )
end

struct NoncollinearPhaseMatch{CT<:NonlinearCrystal} <: PhaseMatch
    refr_data::PMRefractionData
    pm_data::PMNoncollinearData{CT}
    eff_data::PMEfficiencyData
    bw_data::PMBandwidthData
end

function Base.getproperty(ncpm::NoncollinearPhaseMatch, sym::Symbol)
    # Forward PMNoncollinearData and PMRefractionData field requests
    if sym in fieldnames(PMNoncollinearData)
        return getfield(ncpm.pm_data, sym)
    elseif sym in fieldnames(PMRefractionData)
        return getfield(ncpm.refr_data, sym)
    else # Fallback to real fields
        return getfield(ncpm, sym)
    end
end

function NoncollinearPhaseMatch(
    cr::NonlinearCrystal,
    lambda_rrb::NTuple{3,Length},
    temp::Temperature,
    hi_or_lo_rrb::NTuple{3,Symbol},
    theta_pm_rrb::NTuple{3,Angle},
    phi_pm_rrb::NTuple{3,Angle};
    angle_tol::Angle=0.1u"°",
)
    # Sort red lambdas: r1 by definition always has a higher (or equal) wavelength than r2 
    if lambda_rrb[1] < lambda_rrb[2]
        lambda_rrb = lambda_rrb[[2, 1, 3]]
        hi_or_lo_rrb = hi_or_lo_rrb[[2, 1, 3]]
        theta_pm_rrb = theta_pm_rrb[[2, 1, 3]]
        phi_pm_rrb = phi_pm_rrb[[2, 1, 3]]
    end

    rd_rrb = Tuple((RefractionData(hi_or_lo_rrb[i], theta_pm_rrb[i], phi_pm_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(lambda_rrb)))
    refr_data = PMRefractionData(rd_rrb...)

    pm_data = PMNoncollinearData(refr_data, hi_or_lo_rrb, theta_pm_rrb, phi_pm_rrb, cr, lambda_rrb, temp; angle_tol)

    eff_data = PMEfficiencyData(pm_data, refr_data)
    bw_data = PMBandwidthData(pm_data, refr_data)

    return NoncollinearPhaseMatch{typeof(cr)}(refr_data, pm_data, eff_data, bw_data)
end


function Base.show(io::IO, ncpm::NoncollinearPhaseMatch)
    digits = 3

    # Helper
    vec_str(v) = "[" * join(round.(v; digits), ", ") * "]"

    # Header
    @printf(io, "%-29s %s\n", "Crystal:", ncpm.cr.metadata[:description])
    @printf(io, "%-29s %3.2f K (%3.2f °C)\n", "Temperature:",
        ustrip(u"K", ncpm.temp), ustrip(u"°C", ncpm.temp))

    println(io, "────────────────────────────────────────────────────────────────────────────────────────────────────────")

    # Wavelengths
    λs = ustrip.(u"nm", round.(u"nm", ncpm.lambda_rrb; digits))
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Wavelength (nm):", λs...)

    # Polarization types
    if isa(ncpm.cr, UnidirectionalCrystal)
        types = ncpm.pm_data.pm_type
        polar_str = ["$(ncpm.hi_or_lo_rrb[i]) ($(types[1].o_or_e_rrb[i]))" for i in 1:3]
        @printf(io, "%-29s %-25s %-25s %-25s\n", "Type $(types[1].type) PM:", polar_str...)
    else
        @printf(io, "%-29s %-25s %-25s %-25s\n", "Refractive index type:", ncpm.hi_or_lo_rrb...)
        for t in ncpm.pm_data.pm_type
            isnothing(t) && continue
            pols = ["$(t.o_or_e_rrb[i])" for i in 1:3]
            @printf(io, "%-29s %-25s %-25s %-25s\n", "Type $(t.type) PM in $(t.principal_plane) plane:", pols...)
        end
    end

    # Index summary
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Refractive index:",
        auto_fmt.(ncpm.n_rrb; digits)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Group index:",
        auto_fmt.(ncpm.group_index_rrb; digits)...)

    # Walkoff
    w = auto_fmt.(ustrip.(u"mrad", ncpm.walkoff_angle_rrb); digits)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Walkoff angle (mrad):", w...)

    # Directions

    # k angles in degrees
    k_angles_str = [
        "θ: $(round(ustrip(u"°", ncpm.theta_pm_rrb[i]); digits=2))°, ϕ: $(round(ustrip(u"°", ncpm.phi_pm_rrb[i]); digits=2))°"
        for i in 1:3
    ]
    @printf(io, "%-29s %-25s %-25s %-25s\n", "k angle (deg):", k_angles_str...)

    # k directions
    k_dirs = [
        vec_str(angles_to_vector(ncpm.theta_pm_rrb[i], ncpm.phi_pm_rrb[i]))
        for i in 1:3
    ]
    @printf(io, "%-29s %-25s %-25s %-25s\n", "k direction:", k_dirs...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "S direction:", vec_str.(ncpm.S_dir_rrb)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "E direction:", vec_str.(ncpm.E_dir_rrb)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "D direction:", vec_str.(ncpm.D_dir_rrb)...)

    # Dispersion
    gvd = auto_fmt.(ustrip.(u"fs^2/mm", ncpm.beta2_rrb); digits)
    tod = auto_fmt.(ustrip.(u"fs^3/mm", ncpm.beta3_rrb); digits)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "β₂ (fs²/mm):", gvd...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "β₃ (fs³/mm):", tod...)

    # Bandwidths
    ωbw = auto_fmt.(ustrip.(u"GHz*cm", round.(u"GHz*cm", ncpm.bw_data.omega_L_bw; digits)))
    @printf(io, "%-29s %-25s %-25s %-25s\n", "ω BW × L (GHz·cm):", ωbw...)

    @printf(io, "%-29s %-25s\n", "T BW × L (K·cm):",
        auto_fmt(ustrip(u"K*cm", round(u"K*cm", ncpm.bw_data.temp_L_bw; digits))))

    @printf(io, "%-29s %-25s\n", "θ BW × L (mrad·cm):",
        auto_fmt(ustrip(u"mrad*cm", round(u"mrad*cm", ncpm.bw_data.theta_L_bw; digits))))

    @printf(io, "%-29s %-25s\n", "ϕ BW × L (mrad·cm):",
        auto_fmt(ustrip(u"mrad*cm", round(u"mrad*cm", ncpm.bw_data.phi_L_bw; digits))))

    # Efficiency
    @printf(io, "%-29s %-25s\n",
        "d_eff (pm/V):",
        auto_fmt(ustrip(u"pm/V", round(u"pm/V", ncpm.eff_data.d_eff; digits))) *
        " (w/o Miller scaling: " *
        auto_fmt(ustrip(u"pm/V", round(u"pm/V", ncpm.eff_data.d_eff_no_miller; digits))) *
        ")"
    )

    @printf(io, "%-29s %-25s\n",
        "S₀ × L² (W):",
        auto_fmt(ustrip(u"W", round(u"W", ncpm.eff_data.S0_Lsquared; sigdigits=digits))) *
        " (w/o Miller scaling: " *
        auto_fmt(ustrip(u"W", round(u"W", ncpm.eff_data.S0_Lsquared_no_miller; sigdigits=digits))) *
        ")"
    )

    println(io, "────────────────────────────────────────────────────────────────────────────────────────────────────────")
end

function to_yaml(io::IO, ncpm::NoncollinearPhaseMatch)
    # Scalars and vectors (converted to lists to avoid tuple serialization)
    temperature = ustrip(u"K", ncpm.temp)
    wavelength = collect(ustrip.(u"m", ncpm.lambda_rrb))
    refractive_index = collect(ncpm.n_rrb)
    group_index = collect(ncpm.group_index_rrb)
    walkoff_angle = collect(ustrip.(u"rad", ncpm.walkoff_angle_rrb))
    beta2 = collect(ustrip.(u"s^2/m", ncpm.beta2_rrb))
    beta3 = collect(ustrip.(u"s^3/m", ncpm.beta3_rrb))

    # Geometry
    k_angles = [
        Dict(
            "theta" => ustrip(u"rad", ncpm.theta_pm_rrb[i]),
            "phi" => ustrip(u"rad", ncpm.phi_pm_rrb[i])
        ) for i in 1:3
    ]

    k_directions = [
        collect(angles_to_vector(ncpm.theta_pm_rrb[i], ncpm.phi_pm_rrb[i]))
        for i in 1:3
    ]

    geometry = Dict(
        "k_angle" => k_angles,
        "k_direction" => k_directions,
        "S_direction" => [collect(ncpm.S_dir_rrb[i]) for i in 1:3],
        "E_direction" => [collect(ncpm.E_dir_rrb[i]) for i in 1:3],
        "D_direction" => [collect(ncpm.D_dir_rrb[i]) for i in 1:3]
    )

    # Bandwidths
    bandwidths = Dict(
        "omega_L" => collect(ustrip.(u"Hz*m", ncpm.bw_data.omega_L_bw)),
        "temp_L" => ustrip(u"K*m", ncpm.bw_data.temp_L_bw),
        "theta_L" => ustrip(u"rad*m", ncpm.bw_data.theta_L_bw),
        "phi_L" => ustrip(u"rad*m", ncpm.bw_data.phi_L_bw)
    )

    # Efficiency
    efficiency = Dict(
        "d_eff" => ustrip(u"m/V", ncpm.eff_data.d_eff),
        "d_eff_no_miller" => ustrip(u"m/V", ncpm.eff_data.d_eff_no_miller),
        "S0_L2" => ustrip(u"W", ncpm.eff_data.S0_Lsquared),
        "S0_L2_no_miller" => ustrip(u"W", ncpm.eff_data.S0_Lsquared_no_miller)
    )

    # Taylor expansion
    taylor_tpo = [taylor_theta_phi_omega(hi_or_lo, theta_ncpm, phi_ncpm, ncpm.cr, lambda, ncpm.temp) for (hi_or_lo, theta_ncpm, phi_ncpm, lambda) in zip(ncpm.hi_or_lo_rrb, ncpm.theta_pm_rrb, ncpm.phi_pm_rrb, ncpm.lambda_rrb)]
    center_rrb = [t[1] for t in taylor_tpo]
    grad_rrb = [t[2] for t in taylor_tpo]
    hess_rrb = [t[3] for t in taylor_tpo]

    taylor_deriv_theta_phi_omega = Dict(
        "center_rrb" => center_rrb,
        "grad_rrb" => grad_rrb,
        "hess_rrb" => hess_rrb,
    )

    # Top-level dictionary
    data = Dict(
        "crystal" => ncpm.cr.metadata[:description],
        "temperature" => temperature,
        "wavelength" => wavelength,
        "refractive_index" => refractive_index,
        "group_index" => group_index,
        "walkoff_angle" => walkoff_angle,
        "geometry" => geometry,
        "beta2" => beta2,
        "beta3" => beta3,
        "bandwidths" => bandwidths,
        "efficiency" => efficiency,
        "n_taylor_theta_phi_omega" => taylor_deriv_theta_phi_omega
    )

    # Polarization
    if isa(ncpm.cr, UnidirectionalCrystal)
        types = ncpm.pm_data.pm_type
        data["polarization"] = Dict(
            "type" => string(types[1].type),
            "components" => [
                Dict("hi_lo" => ncpm.hi_or_lo_rrb[i], "pol" => types[1].o_or_e_rrb[i])
                for i in 1:3
            ]
        )
    else
        data["refractive_index_type"] = ncpm.hi_or_lo_rrb
        pols = []
        for t in ncpm.pm_data.pm_type
            isnothing(t) && continue
            push!(pols, Dict(
                "type" => string(t.type),
                "principal_plane" => string(t.principal_plane),
                "pols" => t.o_or_e_rrb
            ))
        end
        data["polarization"] = pols
    end

    YAML.write(io, data)
end

"""
    external_ray_dirs(theta_cut, phi_cut, ncpm; n_external_rrb=(1.0, 1.0, 1.0))

Given the cut orientation (`theta_cut`, `phi_cut`) of a birefringent crystal and a noncollinear phasematching configuration `ncpm`, compute the external propagation directions of the three interacting waves after refraction at the crystal surface.

Refraction is computed using Snell's law for each wave individually, taking into account its internal refractive index (`n_internal_rrb`) and a user-specified external refractive index tuple (`n_external_rrb`), e.g., `(1.0, 1.0, 1.0)` for air.
Returns a tuple of three unit vectors representing the external directions of the waves. If total internal reflection occurs for any wave, `nothing` is returned in its place.
"""

function external_ray_dirs(
    theta_cut::Angle,
    phi_cut::Angle,
    ncpm::NoncollinearPhaseMatch;
    n_external_rrb::NTuple{3,<:Real}=(1.0, 1.0, 1.0)
)
    k_dir_internal_rrb = angles_to_vector.(ncpm.theta_pm_rrb, ncpm.phi_pm_rrb)
    n_dir_cut = angles_to_vector(theta_cut, phi_cut)
    n_internal_rrb = ncpm.refr_data.n_rrb

    k_dir_external_rrb = map(zip(k_dir_internal_rrb, n_internal_rrb, n_external_rrb)) do (k_dir_internal, n_i, n_e)
        # Tangential component
        k_tan = k_dir_internal - dot(k_dir_internal, n_dir_cut) * n_dir_cut
        sin_theta_i = norm(k_tan)

        # Snell's law
        sin_theta_e = (n_i / n_e) * sin_theta_i

        if sin_theta_e > 1
            return nothing  # Total internal reflection
        end

        k_tan_ext = (n_i / n_e) * k_tan
        k_n_ext = sqrt(1 - norm(k_tan_ext)^2)
        k_ext = normalize(k_tan_ext + k_n_ext * n_dir_cut)
        return k_ext
    end

    return Tuple(k_dir_external_rrb)
end


function delta_k_noncollinear(
    theta_pm_rrb::NTuple{3,Angle},
    phi_pm_rrb::NTuple{3,Angle},
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    calc_refraction_data::Bool=false,
)
    lambda_rrb = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_rrb])

    if calc_refraction_data
        rd_rrb = Tuple((RefractionData(hi_or_lo_rrb[i], theta_pm_rrb[i], phi_pm_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(lambda_rrb)))
        k_vec_rrb = Tuple((angles_to_vector(rd.theta, rd.phi) * rd.n / rd.lambda for rd in rd_rrb))
    else
        n_rrb = [calc_n_hi_lo(theta_pm_rrb[i], phi_pm_rrb[i], cr, lambda_rrb[i]; temp)[hi_or_lo_rrb[i] === :hi ? 1 : 2] for i in eachindex(lambda_rrb)]
        k_vec_rrb = [angles_to_vector(theta_pm_rrb[i], phi_pm_rrb[i]) * n_rrb[i] / lambda_rrb[i] for i in eachindex(lambda_rrb)]
    end

    delta_k_vec = k_vec_rrb[1] + k_vec_rrb[2] - k_vec_rrb[3]

    return calc_refraction_data ? (delta_k_vec, rd_rrb) : delta_k_vec
end



"""
    delta_k_noncollinear_with_shifting(theta_pm_rrb, phi_pm_rrb, hi_or_lo_rrb, cr; ...)

Same as [`delta_k_noncollinear`](@ref), but allows small parameter perturbations (Δθ, Δϕ, ΔT, Δω_r1, Δω_r2, Δω_b).
The frequency shifts adjust the wavelengths while conserving energy (Δω_r1 + Δω_r2 = Δω_b), enabling bandwidth estimation.

Used internally for computing how Δk varies with perturbations for calculating ΔX · L quantities in [`PMBandwidthData`](@ref).
"""
function delta_k_noncollinear_with_shifting(
    theta_pm_rrb::NTuple{3,Angle},
    phi_pm_rrb::NTuple{3,Angle},
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    delta_theta=0.0u"°",
    delta_phi=0.0u"°",
    delta_temp::Temperature=0.0u"K",
    delta_omega_r1::Frequency=0.0u"GHz",
    delta_omega_r2::Frequency=0.0u"GHz",
    delta_omega_b::Frequency=0.0u"GHz",
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_rrb])

    @assert delta_omega_r1 + delta_omega_r2 == delta_omega_b "Phasematching is only fulfilled for Δω_r1 + Δω_r2 = Δω_b"
    lambda_r1_shifted = shift_lambda_with_freq(lambda_r1, delta_omega_r1)
    lambda_r2_shifted = shift_lambda_with_freq(lambda_r2, delta_omega_r2)
    lambda_b_shifted = shift_lambda_with_freq(lambda_b, delta_omega_b)

    delta_k_noncollinear(
        (theta_pm_rrb .+ delta_theta) .|> u"rad",
        (phi_pm_rrb .+ delta_phi) .|> u"rad",
        hi_or_lo_rrb,
        cr;
        lambda_r1=lambda_r1_shifted,
        lambda_r2=lambda_r2_shifted,
        lambda_b=lambda_b_shifted,
        temp=temp + delta_temp,
    )
end


function temperature_L_bandwidth(pm_data::PMNoncollinearData)
    fun = ΔT -> ustrip.(
        u"m^-1",
        delta_k_noncollinear_with_shifting(
            pm_data.theta_pm_rrb,
            pm_data.phi_pm_rrb,
            pm_data.hi_or_lo_rrb,
            pm_data.cr;
            lambda_r1=pm_data.lambda_rrb[1],
            lambda_r2=pm_data.lambda_rrb[2],
            lambda_b=pm_data.lambda_rrb[3],
            pm_data.temp,
            delta_temp=ΔT * u"K",
        )
    )
    ddeltak_dT = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * K^-1"

    # ΔT_bw * L = 2π / |∂Δk/∂T|
    return 2π / norm(ddeltak_dT) # TODO: Using the norm here is a conservative assumption. A projection in a certain direction could be more suitable here. 
end

function theta_L_bandwidth(pm_data::PMNoncollinearData)
    fun = Δθ -> ustrip.(
        u"m^-1",
        delta_k_noncollinear_with_shifting(
            pm_data.theta_pm_rrb,
            pm_data.phi_pm_rrb,
            pm_data.hi_or_lo_rrb,
            pm_data.cr;
            lambda_r1=pm_data.lambda_rrb[1],
            lambda_r2=pm_data.lambda_rrb[2],
            lambda_b=pm_data.lambda_rrb[3],
            pm_data.temp,
            delta_theta=Δθ * u"rad",
        )
    )
    ddeltak_dtheta = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * rad^-1"

    # Δθ_bw * L = 2π / |∂Δk/∂θ|
    return 2π / norm(ddeltak_dtheta) # TODO: Using the norm here is a conservative assumption. A projection in a certain direction could be more suitable here. 
end

function phi_L_bandwidth(pm_data::PMNoncollinearData)
    fun = Δϕ -> ustrip.(
        u"m^-1",
        delta_k_noncollinear_with_shifting(
            pm_data.theta_pm_rrb,
            pm_data.phi_pm_rrb,
            pm_data.hi_or_lo_rrb,
            pm_data.cr;
            lambda_r1=pm_data.lambda_rrb[1],
            lambda_r2=pm_data.lambda_rrb[2],
            lambda_b=pm_data.lambda_rrb[3],
            pm_data.temp,
            delta_phi=Δϕ * u"rad",
        )
    )
    ddeltak_dphi = ForwardDiff.derivative(
        fun,
        0.0
    ) * 1u"m^-1 * rad^-1"

    # Δϕ_bw * L = 2π / |∂Δk/∂ϕ|
    return 2π / norm(ddeltak_dphi) # TODO: Using the norm here is a conservative assumption. A projection in a certain direction could be more suitable here. 
end