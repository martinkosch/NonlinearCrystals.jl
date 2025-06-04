export PhaseMatch, CollinearPhaseMatch, pm_wavelengths, find_all_pms_along_dimension, find_all_ncpm_over_temp, find_all_ncpm_over_lambda, find_nearest_pm_along_lambda_r_b, find_nearest_pm_along_theta_phi, delta_k

"""
    PMRefractionData

Stores all refractive, group, and dispersion properties for the three interacting waves in a phase-matching configuration.
This includes refractive indices, group indices, walkoff angles, polarization vectors (`E`, `D`, `S`), and dispersion parameters (β₀ to β₃) for the two red and one blue waves (r₁, r₂, b).
This is typically constructed from three [`RefractionData`](@ref) instances.
"""

struct PMRefractionData
    n_rrb::NTuple{3,Float64}
    group_index_rrb::NTuple{3,Float64}
    walkoff_angle_rrb::NTuple{3,typeof(1.0u"rad")}
    D_dir_rrb::NTuple{3,SVector{3,Float64}}
    E_dir_rrb::NTuple{3,SVector{3,Float64}}
    S_dir_rrb::NTuple{3,SVector{3,Float64}}
    beta0_rrb::NTuple{3,typeof(1.0u"m^-1")}
    beta1_rrb::NTuple{3,typeof(1.0u"s * m^-1")}
    beta2_rrb::NTuple{3,typeof(1.0u"s^2 * m^-1")}
    beta3_rrb::NTuple{3,typeof(1.0u"s^3 * m^-1")}
end

function PMRefractionData(rd_rrb::NTuple{3,<:RefractionData})
    # Flip the signs of all field direction vectors (until here chosen randomly) to avoid sign jumps in d_eff
    unified_dir_signs = calc_unified_dir_signs(rd_rrb)

    return PMRefractionData(
        Tuple(r.n for r in rd_rrb),
        Tuple(r.group_index for r in rd_rrb),
        Tuple(r.walkoff_angle for r in rd_rrb),
        Tuple(r.D_dir for r in rd_rrb) .* unified_dir_signs,
        Tuple(r.E_dir for r in rd_rrb) .* unified_dir_signs,
        Tuple(r.S_dir for r in rd_rrb),
        Tuple(r.beta0 for r in rd_rrb),
        Tuple(r.beta1 for r in rd_rrb),
        Tuple(r.beta2 for r in rd_rrb),
        Tuple(r.beta3 for r in rd_rrb),
    )
end

PMRefractionData(rd_r1::RefractionData, rd_r2::RefractionData, rd_b::RefractionData) = PMRefractionData((rd_r1, rd_r2, rd_b))

"""
    PMType

Describes the polarization type of a phase-matching configuration, including the principal plane (e.g. `:XY` or `:UD` for uniaxial crystals),
the polarization of each wave (`:o` or `:e`), and a human-readable type classification like `"I"`, `"II/III"`, `"IV"`, or `"V"`.

Used for diagnostic purposes and auto-labeling of phase-match types.
"""
struct PMType
    principal_plane::Symbol
    o_or_e_rrb::NTuple{3,Symbol}
    type::String
end

function Base.show(io::IO, pmt::PMType)
    if pmt.principal_plane === :UD
        @printf(io, "%-32s %-3s %-3s %-3s\n", "Unidirectional type $(pmt.type) PM:", pmt.o_or_e_rrb...)
    else
        pols = ["$(pmt.o_or_e_rrb[i])" for i in 1:3]
        @printf(io, "%-32s %-3s %-3s %-3s\n", "Bidirectional type $(pmt.type) PM in $(pmt.principal_plane) plane:", pols...)
    end
end

function PMType(
    principal_plane::Symbol,
    refr_data::PMRefractionData;
)
    o_or_e_rrb = Tuple(assign_o_or_e(principal_plane, E_dir) for E_dir in refr_data.E_dir_rrb)

    if o_or_e_rrb[3] === :e
        if o_or_e_rrb[1] === o_or_e_rrb[2] === :e # :e, :e, :e
            type = "IV"
        elseif o_or_e_rrb[1] === o_or_e_rrb[2] === :o # :o, :o, :e
            type = "I"
        else # :o, :e, :e or e:, :o, :e
            type = "II/III"
        end
    else # o_or_e_rrb[3] === :o
        if o_or_e_rrb[1] === o_or_e_rrb[2] === :o # :o, :o, :o
            type = "V"
        elseif o_or_e_rrb[1] === o_or_e_rrb[2] === :e # :e, :e, :o
            type = "I" # Extended scheme: "I/VIII"
        else # :o, :e, :o or e:, :o, :o
            type = "II/III" # Extended scheme: "VI/VII"
        end
    end

    return PMType(principal_plane, o_or_e_rrb, type)
end

"""
    PMCollinearData

Represents the geometry and wavelengths of a collinear phase-matching configuration, including propagation angles `theta_pm` and `phi_pm`,
crystal reference, temperature, and polarization type data (`hi_or_lo_rrb` and optional [`PMType`](@ref) classification).

This is the geometric input used to calculate [`CollinearPhaseMatch`](@ref) properties.
"""
struct PMCollinearData{CT<:NonlinearCrystal}
    hi_or_lo_rrb::NTuple{3,Symbol}
    pm_type::NTuple{2,Union{Nothing,PMType}}
    theta_pm::typeof(1.0u"rad")
    phi_pm::typeof(1.0u"rad")
    cr::CT
    lambda_rrb::NTuple{3,typeof(1.0u"m")}
    temp::typeof(1.0u"K")
end

function PMCollinearData(
    refr_data::PMRefractionData,
    hi_or_lo_rrb::NTuple{3,Symbol},
    theta_pm::Angle,
    phi_pm::Angle,
    cr::NonlinearCrystal,
    lambda_rrb::NTuple{3,Length},
    temp::Temperature;
    angle_tol::Angle=0.1u"°",
)
    if isa(cr, UnidirectionalCrystal)
        # o and e are valid classifiers for UnidirectionalCrystal even if propagating outside of the pricipal planes
        pm_type = (PMType(:UD, refr_data), nothing)
    else
        principal_planes = find_principal_planes(theta_pm, phi_pm; angle_tol)
        pm_type = Tuple((isnothing(p) ? nothing : PMType(p, refr_data) for p in principal_planes))
    end

    return PMCollinearData{typeof(cr)}(
        hi_or_lo_rrb,
        pm_type,
        theta_pm,
        phi_pm,
        cr,
        lambda_rrb,
        temp,
    )
end

"""
    PMEfficiencyData

Stores the effective nonlinear coefficient `d_eff` (with and without Miller scaling) and the associated conversion factor `S₀·L²`.
Used to compare the efficiency of different phase-matching solutions based on geometry and polarization.
"""
struct PMEfficiencyData
    d_eff::typeof(1.0u"m/V")
    d_eff_no_miller::typeof(1.0u"m/V")
    S0_Lsquared::typeof(1.0u"W")
    S0_Lsquared_no_miller::typeof(1.0u"W")
end

function PMEfficiencyData(pm_data::PMCollinearData, refr_data::PMRefractionData)
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

"""
    PMBandwidthData

Represents the phase-matching bandwidths of a collinear configuration, defined as the range of parameters over which the phase mismatch Δk remains within ±π / L for a given interaction length L.

The stored values are ΔX · L quantities (units of Hz·m, K·m, rad·m), which characterize the *tolerance* of the phase-matching condition to changes in each parameter X:

- `omega_L_bw`: Tuple of angular frequency bandwidths (Δω · L) for the three waves (r₁, r₂, b). Each entry is computed by holding that wave's frequency fixed while adjusting the other two to maintain the phase-matching condition (1/λ_r₁ + 1/λ_r₂ = 1/λ_b).  
- `temp_L_bw`: Temperature bandwidth (ΔT · L) computed by evaluating how Δk varies with small temperature changes, with all wavelengths and geometry fixed.
- `theta_L_bw`: Angular bandwidth with respect to polar angle θ (Δθ · L), evaluated at fixed temperature and wavelengths.
- `phi_L_bw`: Angular bandwidth with respect to azimuthal angle ϕ (Δϕ · L), also at fixed conditions.

These values describe how sensitive the phase-matching is to deviations in each parameter. For example, a small `theta_L_bw` means tight angular alignment is needed, while a large `temp_L_bw` suggests thermally robust operation.

Bandwidths are calculated as:
```math
ΔX · L = \\frac{2π}{|∂Δk/∂X|}
```
where Δk is the phase mismatch, and X is one of frequency, temperature, or angle. 
"""
struct PMBandwidthData
    omega_L_bw::NTuple{3,typeof(1.0u"Hz*m")}
    temp_L_bw::typeof(1.0u"K*m")
    theta_L_bw::typeof(1.0u"rad*m")
    phi_L_bw::typeof(1.0u"rad*m")
end

function PMBandwidthData(pm_data::PMCollinearData, refr_data::PMRefractionData)
    omega_L_bw = lambda_L_bandwidths(refr_data)
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

abstract type PhaseMatch end

"""
    CollinearPhaseMatch <: PhaseMatch

Encapsulates a full collinear phase-matching solution, including refractive data, geometric parameters, 
nonlinear efficiency, and bandwidths (fields `refr_data`, `pm_data`, `eff_data`, and `bw_data`, respectively).
In addition, field access is forwarded to the internal fields of both `refr_data` and `pm_data` for convenience.
"""
struct CollinearPhaseMatch{CT<:NonlinearCrystal} <: PhaseMatch
    refr_data::PMRefractionData
    pm_data::PMCollinearData{CT}
    eff_data::PMEfficiencyData
    bw_data::PMBandwidthData
end

function Base.getproperty(cpm::CollinearPhaseMatch, sym::Symbol)
    # Forward PMCollinearData and PMRefractionData field requests
    if sym in fieldnames(PMCollinearData)
        return getfield(cpm.pm_data, sym)
    elseif sym in fieldnames(PMRefractionData)
        return getfield(cpm.refr_data, sym)
    else # Fallback to real fields
        return getfield(cpm, sym)
    end
end

function CollinearPhaseMatch(
    cr::NonlinearCrystal,
    lambda_rrb::NTuple{3,Length},
    temp::Temperature,
    hi_or_lo_rrb::NTuple{3,Symbol},
    theta_pm::Angle,
    phi_pm::Angle;
    angle_tol::Angle=0.1u"°",
)
    # Sort red lambdas: r1 by definition always has a higher (or equal) wavelength than r2 
    if lambda_rrb[1] < lambda_rrb[2]
        lambda_rrb = lambda_rrb[[2, 1, 3]]
        hi_or_lo_rrb = hi_or_lo_rrb[[2, 1, 3]]
    end

    rd_rrb = [RefractionData(hi_or_lo_rrb[i], theta_pm, phi_pm, cr, lambda_rrb[i]; temp) for i in eachindex(lambda_rrb)]
    refr_data = PMRefractionData(rd_rrb...)

    pm_data = PMCollinearData(refr_data, hi_or_lo_rrb, theta_pm, phi_pm, cr, lambda_rrb, temp; angle_tol)

    eff_data = PMEfficiencyData(pm_data, refr_data)
    bw_data = PMBandwidthData(pm_data, refr_data)

    return CollinearPhaseMatch{typeof(cr)}(refr_data, pm_data, eff_data, bw_data)
end

function calc_unified_dir_signs(
    rd_rrb::NTuple{3,<:RefractionData}
)
    E_dir_rrb = Tuple(r.E_dir for r in rd_rrb)
    first_optical_axis_b = angles_to_vector(optical_axis_angle(rd_rrb[3].cr, rd_rrb[3].lambda, rd_rrb[3].temp), 0.0u"°")

    # TODO: Check this and the signs of d_eff at opposite hemispheres
    # This sign correction is supposed to prevent sign jumps, e.g. in d_eff
    d_b_oa = dot(E_dir_rrb[3], first_optical_axis_b)
    d_b_ref = dot(E_dir_rrb[3], cross(rd_rrb[3].S_dir, first_optical_axis_b))
    sign_b = (abs(d_b_oa) > abs(d_b_ref)) ? sign(d_b_oa) : sign(d_b_ref)

    dir_ref = cross(E_dir_rrb[3], rd_rrb[3].S_dir)
    d_r1_ref = dot(dir_ref, E_dir_rrb[1])
    d_r1_b = dot(E_dir_rrb[3], E_dir_rrb[1])
    d_r2_ref = dot(dir_ref, E_dir_rrb[2])
    d_r2_b = dot(E_dir_rrb[3], E_dir_rrb[2])

    sign_r1 = (abs(d_r1_ref) > abs(d_r1_b)) ? sign(d_r1_ref) : sign(d_r1_b)
    sign_r2 = (abs(d_r2_ref) > abs(d_r2_b)) ? sign(d_r2_ref) : sign(d_r2_b)
    return (sign_r1, sign_r2, sign_b)
end

function Base.show(io::IO, cpm::CollinearPhaseMatch)
    digits = 3

    # Helper
    vec_str(v) = "[" * join(round.(v; digits), ", ") * "]"

    # Header
    @printf(io, "%-29s %s\n", "Crystal:", cpm.cr.metadata[:description])
    @printf(io, "%-29s θ: %3.2f°, ϕ: %3.2f°\n", "k angles:",
        ustrip(u"°", cpm.theta_pm), ustrip(u"°", cpm.phi_pm))
    @printf(io, "%-29s %-25s\n", "k direction:", vec_str(angles_to_vector(cpm.theta_pm, cpm.phi_pm)))
    @printf(io, "%-29s %3.2f K (%3.2f °C)\n", "Temperature:",
        ustrip(u"K", cpm.temp), ustrip(u"°C", cpm.temp))

    println(io, "────────────────────────────────────────────────────────────────────────────────────────────────────────")

    # Wavelengths
    λs = ustrip.(u"nm", round.(u"nm", cpm.lambda_rrb; digits))
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Wavelength (nm):", λs...)

    # Polarization types
    if isa(cpm.cr, UnidirectionalCrystal)
        types = cpm.pm_data.pm_type
        polar_str = ["$(cpm.hi_or_lo_rrb[i]) ($(types[1].o_or_e_rrb[i]))" for i in 1:3]
        @printf(io, "%-29s %-25s %-25s %-25s\n", "Type $(types[1].type) PM:", polar_str...)
    else
        @printf(io, "%-29s %-25s %-25s %-25s\n", "Refractive index type:", cpm.hi_or_lo_rrb...)
        for t in cpm.pm_data.pm_type
            isnothing(t) && continue
            pols = ["$(t.o_or_e_rrb[i])" for i in 1:3]
            @printf(io, "%-29s %-25s %-25s %-25s\n", "Type $(t.type) PM in $(t.principal_plane) plane:", pols...)
        end
    end

    # Index summary
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Refractive index:",
        auto_fmt.(cpm.n_rrb; digits)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Group index:",
        auto_fmt.(cpm.group_index_rrb; digits)...)

    # Walkoff
    w = auto_fmt.(ustrip.(u"mrad", cpm.walkoff_angle_rrb); digits)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "Walkoff angle (mrad):", w...)

    # Directions
    @printf(io, "%-29s %-25s %-25s %-25s\n", "S direction:", vec_str.(cpm.S_dir_rrb)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "E direction:", vec_str.(cpm.E_dir_rrb)...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "D direction:", vec_str.(cpm.D_dir_rrb)...)

    # Dispersion
    gvd = auto_fmt.(ustrip.(u"fs^2/mm", cpm.beta2_rrb); digits)
    tod = auto_fmt.(ustrip.(u"fs^3/mm", cpm.beta3_rrb); digits)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "β₂ (fs²/mm):", gvd...)
    @printf(io, "%-29s %-25s %-25s %-25s\n", "β₃ (fs³/mm):", tod...)

    # Bandwidths
    ωbw = auto_fmt.(ustrip.(u"GHz*cm", round.(u"GHz*cm", cpm.bw_data.omega_L_bw; digits)))
    @printf(io, "%-29s %-25s %-25s %-25s\n", "ω BW × L (GHz·cm):", ωbw...)

    @printf(io, "%-29s %-25s\n", "T BW × L (K·cm):",
        auto_fmt(ustrip(u"K*cm", round(u"K*cm", cpm.bw_data.temp_L_bw; digits))))

    @printf(io, "%-29s %-25s\n", "θ BW × L (mrad·cm):",
        auto_fmt(ustrip(u"mrad*cm", round(u"mrad*cm", cpm.bw_data.theta_L_bw; digits))))

    @printf(io, "%-29s %-25s\n", "ϕ BW × L (mrad·cm):",
        auto_fmt(ustrip(u"mrad*cm", round(u"mrad*cm", cpm.bw_data.phi_L_bw; digits))))

    # Efficiency
    @printf(io, "%-29s %-25s\n",
        "d_eff (pm/V):",
        auto_fmt(ustrip(u"pm/V", round(u"pm/V", cpm.eff_data.d_eff; digits))) *
        " (w/o Miller scaling: " *
        auto_fmt(ustrip(u"pm/V", round(u"pm/V", cpm.eff_data.d_eff_no_miller; digits))) *
        ")"
    )

    @printf(io, "%-29s %-25s\n",
        "S₀ × L² (W):",
        auto_fmt(ustrip(u"W", round(u"W", cpm.eff_data.S0_Lsquared; sigdigits=digits))) *
        " (w/o Miller scaling: " *
        auto_fmt(ustrip(u"W", round(u"W", cpm.eff_data.S0_Lsquared_no_miller; sigdigits=digits))) *
        ")"
    )

    println(io, "────────────────────────────────────────────────────────────────────────────────────────────────────────")
end

"""
    pm_wavelengths(; lambda_r1=nothing, lambda_r2=nothing, lambda_b=nothing)

Given any two of the three wavelengths involved in a second-order three-wave interaction (`λ_r1`, `λ_r2`, `λ_b`), compute the third to satisfy energy conservation:
```math
1 / λ_{b} = 1 / λ_{r,1} + 1 / λ_{r,2}
```

Returns the tuple (λ_r1, λ_r2, λ_b). Throws an error if the relation is not fulfilled or fewer than two inputs are provided. 
"""
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
    theta_pm::Angle,
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::UnidirectionalCrystal;
    kwargs...
)
    # Use symmetry in unidirectional crystals
    return delta_k(theta_pm, 0.0u"°", hi_or_lo_rrb, cr; kwargs...)
end

"""
    delta_k(theta_pm, phi_pm, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp)

Computes the phase mismatch Δk = k_r1 + k_r2 - k_b for a given configuration.
Each `k_i` is calculated from `2π * n_i / λ_i` based on the given crystal `cr`, propagation angles `theta_pm`, `phi_pm`, and the three polarization types `hi_or_lo_rrb`.
"""
function delta_k(
    theta_pm::Angle,
    phi_pm::Angle,
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr)
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    @assert all([p in [:hi, :lo] for p in hi_or_lo_rrb])

    n_r1 = calc_n_hi_lo(theta_pm, phi_pm, cr, lambda_r1; temp)[hi_or_lo_rrb[1] === :hi ? 1 : 2]
    n_r2 = calc_n_hi_lo(theta_pm, phi_pm, cr, lambda_r2; temp)[hi_or_lo_rrb[2] === :hi ? 1 : 2]
    n_b = calc_n_hi_lo(theta_pm, phi_pm, cr, lambda_b; temp)[hi_or_lo_rrb[3] === :hi ? 1 : 2]

    return 2π * (n_r1 / lambda_r1 + n_r2 / lambda_r2 - n_b / lambda_b)
end

"""
    delta_k_with_shifting(theta_pm, phi_pm, hi_or_lo_rrb, cr; ...)

Same as [`delta_k`](@ref), but allows small parameter perturbations (Δθ, Δϕ, ΔT, Δω_r1, Δω_r2, Δω_b).
The frequency shifts adjust the wavelengths while conserving energy (Δω_r1 + Δω_r2 = Δω_b), enabling bandwidth estimation.

Used internally for computing how Δk varies with perturbations for calculating ΔX · L quantities in [`PMBandwidthData`](@ref).
"""
function delta_k_with_shifting(
    theta_pm::Angle,
    phi_pm::Angle,
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

    delta_k(
        (theta_pm + delta_theta) |> u"rad",
        (phi_pm + delta_phi) |> u"rad",
        hi_or_lo_rrb,
        cr;
        lambda_r1=lambda_r1_shifted,
        lambda_r2=lambda_r2_shifted,
        lambda_b=lambda_b_shifted,
        temp=temp + delta_temp,
    )
end

function calc_delta_k_map(
    theta::Angle,
    phi::Angle,
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_b_min::Union{Nothing,Length}=nothing,
    lambda_b_max::Union{Nothing,Length}=nothing,
    lambda_r12_min::Union{Nothing,Length}=nothing,
    lambda_r12_max::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=100,
)
    isnothing(lambda_b_min) && (lambda_b_min = valid_lambda_range(cr)[1])
    isnothing(lambda_b_max) && (lambda_b_max = valid_lambda_range(cr)[2] / 2)
    isnothing(lambda_r12_min) && (lambda_r12_min = valid_lambda_range(cr)[1])
    isnothing(lambda_r12_max) && (lambda_r12_max = valid_lambda_range(cr)[2])
    range_lambda_b = LinRange(lambda_b_min, lambda_b_max, ngrid) .|> u"µm"
    range_lambda_r12 = LinRange(lambda_r12_min, lambda_r12_max, ngrid) .|> u"µm"

    all_delta_k = zeros(length(range_lambda_b), length(range_lambda_r12)) * u"m^-1"
    for i in CartesianIndices(all_delta_k)
        (i_b, i_r12) = (i[1], i[2])
        lambda_b = range_lambda_b[i_b]
        lambda_r12 = range_lambda_r12[i_r12]
        lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_b, lambda_r1=lambda_r12)
        if all(is_lambda_valid.([lambda_r1, lambda_r2, lambda_b], cr))
            all_delta_k[i_b, i_r12] = delta_k(theta, phi, hi_or_lo_rrb, cr; temp, lambda_r1, lambda_r2, lambda_b)
        else
            all_delta_k[i_b, i_r12] = NaN * u"m^-1"
        end
    end
    return range_lambda_r12, range_lambda_b, all_delta_k
end

function lambda_L_bandwidths(refr_data::PMRefractionData)
    # ∂Δk/∂f_r2 = ∂k_r2/∂f_r2 - ∂k_b/∂f_b for f_r1=const.
    ddelta_k_df_r2_neg_b = 2π / c_0 * (refr_data.group_index_rrb[2] - refr_data.group_index_rrb[3])

    # ∂Δk/∂f_r1 = ∂k_r1/∂f_r1 - ∂k_b/∂f_b for f_r2=const.
    ddelta_k_df_r1_neg_b = 2π / c_0 * (refr_data.group_index_rrb[1] - refr_data.group_index_rrb[3])

    # ∂Δk/∂f_r1 = ∂k_r1/∂f_r1 - ∂k_r2/∂f_r2 for f_b=const.
    ddelta_k_df_r1_neg_r2 = 2π / c_0 * (refr_data.group_index_rrb[1] - refr_data.group_index_rrb[2])

    # Δf_k_bw  * L = 2π / |∂k_i/∂f_i - ∂k_j/∂f_j| = 2π / |1 / v_g_i - 1 / v_g_j|
    return (
        2π / abs(ddelta_k_df_r2_neg_b),
        2π / abs(ddelta_k_df_r1_neg_b),
        2π / abs(ddelta_k_df_r1_neg_r2)
    )
end

function temperature_L_bandwidth(pm_data::PMCollinearData)
    fun = ΔT -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            pm_data.theta_pm,
            pm_data.phi_pm,
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
    return 2π / abs(ddeltak_dT)
end

function theta_L_bandwidth(pm_data::PMCollinearData)
    fun = Δθ -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            pm_data.theta_pm,
            pm_data.phi_pm,
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
    return 2π / abs(ddeltak_dtheta)
end

function phi_L_bandwidth(pm_data::PMCollinearData)
    fun = Δϕ -> ustrip(
        u"m^-1",
        delta_k_with_shifting(
            pm_data.theta_pm,
            pm_data.phi_pm,
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
    return 2π / abs(ddeltak_dphi)
end


"""
    find_nearest_pm_along_theta_phi(theta_target, phi_target, pol_rrb, cr; ...)

Searches for the phase-matching configuration closest to a given direction (`theta_target`, `phi_target`) by 
scanning phase-matches at fixed θ or fixed ϕ, using `find_all_pms_along_dimension`.

The search returns the [`CollinearPhaseMatch`](@ref) whose propagation direction is most aligned (in dot product) with the target vector.
Returns `nothing` if no solution is found.
"""
function find_nearest_pm_along_theta_phi(
    theta_target::Angle,
    phi_target::Angle,
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    all_pm_θ_fixed = find_all_pms_along_dimension(pol_rrb, cr; lambda_r1_fixed=lambda_r1, lambda_r2_fixed=lambda_r2, lambda_b_fixed=lambda_b, temp_min=temp, temp_max=temp, ngrid, tol, theta_fixed=theta_target)
    all_pm_ϕ_fixed = find_all_pms_along_dimension(pol_rrb, cr; lambda_r1_fixed=lambda_r1, lambda_r2_fixed=lambda_r2, lambda_b_fixed=lambda_b, temp_min=temp, temp_max=temp, ngrid, tol, phi_fixed=phi_target)

    all_pm_candidates = [[pm for pm in all_pm_θ_fixed]; [pm for pm in all_pm_ϕ_fixed]]
    isempty(all_pm_candidates) && return nothing

    all_theta_pm_candidates = [pm.theta_pm for pm in all_pm_candidates]
    all_phi_pm_candidates = [pm.phi_pm for pm in all_pm_candidates]

    target_vec = angles_to_vector(theta_target, phi_target)
    pm_candidate_vecs = angles_to_vector.(all_theta_pm_candidates, all_phi_pm_candidates)
    i_nearest = findmax([dot(pm_vec, target_vec) for pm_vec in pm_candidate_vecs])[2]

    return all_pm_candidates[i_nearest]
end

"""
    find_nearest_pm_along_lambda_r_b(pol_rrb, cr; lambda_r1, lambda_r2, lambda_b, ...)

Searches for the best phase-match for a given pair of red and blue wavelengths. You must specify exactly one of `lambda_r1` or `lambda_r2`, plus `lambda_b`.

Scans along both the missing red wavelength and the blue wavelength, finds all matching configurations (via [`find_all_pms_along_dimension`](@ref)),
and returns the one closest to the provided values.

Used to fine-tune a target wavelength combination to an actual valid phase-match geometry.
"""
function find_nearest_pm_along_lambda_r_b(
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp::Temperature=default_temp(cr),
    principal_axis::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    @assert count(isnothing.([lambda_r1, lambda_r2])) == 1 "Either lambda_r1 or lambda_r2 must be given and the other one must be nothing."

    # Search along lambda_b
    all_pm_lambda_b = find_all_pms_along_dimension(pol_rrb, cr;
        lambda_b_fixed=lambda_b,
        temp_min=temp,
        temp_max=temp,
        principal_axis,
        ngrid=ngrid,
        tol=tol
    )

    # Search along lambda_r1 or lambda_r2
    all_pm_lambda_r = find_all_pms_along_dimension(pol_rrb, cr;
        lambda_r1_fixed=lambda_r1,
        lambda_r2_fixed=lambda_r2,
        temp_min=temp,
        temp_max=temp,
        principal_axis,
        ngrid=ngrid,
        tol=tol
    )

    # Combine candidates
    all_pm_candidates = vcat(all_pm_lambda_b, all_pm_lambda_r)
    isempty(all_pm_candidates) && return nothing

    # Compute distance to target wavelength(s)
    distances = map(all_pm_candidates) do pm
        λs = pm.lambda_rrb
        λ_r = isnothing(lambda_r1) ? λs[2] : λs[1]  # whichever was scanned
        λ_b = λs[3]
        Δλ_r = abs(λ_r - (lambda_r1 === nothing ? lambda_r2 : lambda_r1))
        Δλ_b = abs(λ_b - lambda_b)
        Δλ_r^2 + Δλ_b^2
    end

    i_nearest = findmin(distances)[2]
    return all_pm_candidates[i_nearest]
end

"""
    find_all_ncpm_over_temp(pol_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp_min, temp_max, ...)

Returns all noncritical phase-matching (NCPM) configurations where Δk = 0 can be achieved by tuning the *temperature*.
You must specify at least two of the three wavelengths, and a temperature range `[temp_min, temp_max]`. The search is limited to propagation along principal axes.
Returns an array of [`CollinearPhaseMatch`](@ref) solutions. The return value can be an empty array if no solution is found.
"""
function find_all_ncpm_over_temp(
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    temp_min::Temperature=default_temp(cr),
    temp_max::Temperature=600u"K",
    principal_axis::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    return find_all_pms_along_dimension(
        pol_rrb,
        cr;
        lambda_r1_fixed=lambda_r1,
        lambda_r2_fixed=lambda_r2,
        lambda_b_fixed=lambda_b,
        temp_min,
        temp_max,
        principal_axis,
        ngrid,
        tol,
    )
end

"""
    find_all_ncpm_over_lambda(pol_rrb, cr, temp; lambda_r1, lambda_r2, lambda_b, ...)

Finds all noncritical phase-matches (NCPM) by varying wavelength, while keeping the temperature fixed.
You must specify exactly one of the three wavelengths (`lambda_r1`, `lambda_r2`, or `lambda_b`) and leave the others as `nothing`. The scan range is chosen automatically based on the crystal's validity domain.
Returns a list of [`CollinearPhaseMatch`](@ref) solutions found within the allowed range. The return value can be an empty array if no solution is found.
"""
function find_all_ncpm_over_lambda(
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal,
    temp::Temperature=default_temp(cr);
    lambda_r1::Union{Nothing,Length}=nothing,
    lambda_r2::Union{Nothing,Length}=nothing,
    lambda_b::Union{Nothing,Length}=nothing,
    principal_axis::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    @assert count(isnothing.([lambda_r1, lambda_r2, lambda_b])) == 2 "Only one wavelength must be specified."
    return find_all_pms_along_dimension(
        pol_rrb,
        cr;
        lambda_r1_fixed=lambda_r1,
        lambda_r2_fixed=lambda_r2,
        lambda_b_fixed=lambda_b,
        temp_min=temp,
        temp_max=temp,
        principal_axis,
        ngrid,
        tol,
    )
end

"""
    find_all_pms_along_dimension(pol_rrb, cr; ...)

Performs a one-dimensional search to find all collinear phase-matching solutions (Δk ≈ 0) for a given polarization tuple `pol_rrb`.

Depending on which parameters are held fixed (`lambda_r1`, `lambda_r2`, `lambda_b`, `temp_min == temp_max`, or angular constraints like `theta_fixed`), this scans over:

- Wavelength (λ_b or λ_r1)
- Temperature
- Angle (θ or ϕ)

This function is the main backend used by higher-level routines like [`find_nearest_pm_along_theta_phi`](@ref) or [`find_all_ncpm_over_temp`](@ref).

Returns a list of [`CollinearPhaseMatch`](@ref) solutions found by detecting sign changes in Δk via a global search and refining them with bisection.
"""
function find_all_pms_along_dimension(
    pol_rrb::NTuple{3,Symbol},
    cr::NonlinearCrystal;
    lambda_r1_fixed::Union{Nothing,Length}=nothing,
    lambda_r2_fixed::Union{Nothing,Length}=nothing,
    lambda_b_fixed::Union{Nothing,Length}=nothing,
    temp_min::Temperature=default_temp(cr),
    temp_max::Temperature=default_temp(cr),
    principal_axis::Union{AbstractVector{Symbol},Symbol,Nothing}=nothing,
    theta_fixed=nothing,
    phi_fixed=nothing,
    ngrid=500,
    tol=1e-14u"nm^-1",
)
    λmin, λmax = valid_lambda_range(cr)
    temp_fixed = (temp_min == temp_max) ? temp_min : nothing

    if !isnothing(theta_fixed) || !isnothing(phi_fixed)
        @assert isnothing(principal_axis) "principal_axis must be `nothing` if searching along fixed θ or ϕ directions."
        @assert !isnothing(temp_fixed) "Fix temperature by selecting temp_min = temp_max if searching along fixed θ or ϕ directions."
        return _pms_vs_angle(pol_rrb, cr; lambda_r1=lambda_r1_fixed, lambda_r2=lambda_r2_fixed, lambda_b=lambda_b_fixed, temp=temp_fixed, theta_fixed, phi_fixed, ngrid, tol)
    end

    red_fixed_count = count(.!isnothing.([lambda_r1_fixed, lambda_r2_fixed]))
    @assert (count(.!isnothing.([temp_fixed, lambda_b_fixed])) + red_fixed_count) == 2 "You must provide exactly two of those: (temp_min == temp_max), lambda_b_fixed, lambda_r1_fixed, and lambda_r2_fixed."

    hi_or_lo_rrb = pol_rrb # TODO: Enable searching along temperature and wavelengths with :o/:e notation

    if isnothing(principal_axis) && isnothing(theta_fixed) && isnothing(phi_fixed)
        principal_axis = [:X, :Y, :Z]
    end
    θ_ϕ_principal_axes = axes_to_θ_ϕ(principal_axis)

    if isnothing(temp_fixed)
        return _pms_vs_temperature(θ_ϕ_principal_axes, hi_or_lo_rrb, cr, lambda_r1_fixed, lambda_r2_fixed, lambda_b_fixed, temp_min, temp_max, ngrid, tol)
    elseif isnothing(lambda_b_fixed)
        return _pms_vs_lambda_b(θ_ϕ_principal_axes, hi_or_lo_rrb, cr, lambda_r1_fixed, lambda_r2_fixed, temp_fixed, λmin, λmax, ngrid, tol)
    elseif isnothing(lambda_r1_fixed) && isnothing(lambda_r2_fixed)
        return _pms_vs_lambda_r1(θ_ϕ_principal_axes, hi_or_lo_rrb, cr, lambda_b_fixed, temp_fixed, λmin, λmax, ngrid, tol)
    else
        error("This should never be reached.")
    end
end

function _pms_vs_angle(pol_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp, theta_fixed, phi_fixed, ngrid, tol)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)
    lambda_rrb = (lambda_r1, lambda_r2, lambda_b)

    # Transform :o/e into :hi/:lo if possible
    hi_or_lo_rrb = polarization_rrb_to_hi_lo(pol_rrb, cr, lambda_rrb; temp)

    pms = CollinearPhaseMatch[]

    if !isnothing(theta_fixed) && isnothing(phi_fixed)
        all_ϕ = range(0, 2π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(theta_fixed, ϕ, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp) for ϕ in all_ϕ]

        for i in 1:(length(all_ϕ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0
                ϕ_sol = find_zero(ϕ -> delta_k(theta_fixed, ϕ, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_ϕ[i], all_ϕ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_rrb, temp, hi_or_lo_rrb, theta_fixed, ϕ_sol))
            end
        end

    elseif isnothing(theta_fixed) && !isnothing(phi_fixed)
        all_θ = range(0, π, length=ngrid) * u"rad"
        all_delta_k = [delta_k(θ, phi_fixed, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp) for θ in all_θ]

        for i in 1:(length(all_θ)-1)
            if ustrip(all_delta_k[i] * all_delta_k[i+1]) < 0
                θ_sol = find_zero(θ -> delta_k(θ, phi_fixed, hi_or_lo_rrb, cr; lambda_r1, lambda_r2, lambda_b, temp),
                    (all_θ[i], all_θ[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, lambda_rrb, temp, hi_or_lo_rrb, θ_sol, phi_fixed))
            end
        end
    else
        error("You must provide either theta_fixed or phi_fixed (but not both).")
    end

    return pms
end

function _pms_vs_temperature(θ_ϕs, hi_or_lo_rrb, cr, λr1_fix, λr2_fix, λb_fix, temp_min, temp_max, ngrid, tol)
    pms = CollinearPhaseMatch[]
    λr1, λr2, λb = pm_wavelengths(; lambda_r1=λr1_fix, lambda_r2=λr2_fix, lambda_b=λb_fix)

    for θ_ϕ in θ_ϕs
        temps = range(temp_min, temp_max, length=ngrid)
        Δk = [delta_k(θ_ϕ..., hi_or_lo_rrb, cr; temp, lambda_r1=λr1, lambda_r2=λr2, lambda_b=λb) for temp in temps]

        for i in 1:length(temps)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                fun = T -> delta_k(θ_ϕ..., hi_or_lo_rrb, cr; temp=T, lambda_r1=λr1, lambda_r2=λr2, lambda_b=λb)
                T_sol = find_zero(fun, (temps[i], temps[i+1]), Bisection(), atol=tol)
                push!(pms, CollinearPhaseMatch(cr, (λr1, λr2, λb), T_sol, hi_or_lo_rrb, θ_ϕ...))
            end
        end
    end

    return pms
end

function _pms_vs_lambda_b(θ_ϕs, hi_or_lo_rrb, cr, λr1_fix, λr2_fix, temp, λmin, λmax, ngrid, tol)
    pms = CollinearPhaseMatch[]

    for θ_ϕ in θ_ϕs
        if !isnothing(λr1_fix)
            if 1 / (1 / λr1_fix + 1 / λmax) < λmin
                λb_range = []
            else
                λb_range = range(
                    λmin,
                    1 / (1 / λr1_fix + 1 / λmax),
                    length=ngrid
                )
                Δk = [delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r1=λr1_fix, lambda_b=λb, temp=temp) for λb in λb_range]
            end
        else # !isnothing(λr2_fix)
            if 1 / (1 / λr2_fix + 1 / λmax) < λmin
                λb_range = []
            else
                λb_range = range(
                    λmin,
                    1 / (1 / λr2_fix + 1 / λmax),
                    length=ngrid
                )
                Δk = [delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r2=λr2_fix, lambda_b=λb, temp=temp) for λb in λb_range]
            end
        end

        for i in 1:length(λb_range)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                if !isnothing(λr1_fix)
                    fun = λb -> delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r1=λr1_fix, lambda_b=λb, temp=temp)
                    λb_sol = find_zero(fun, (λb_range[i], λb_range[i+1]), Bisection(), atol=tol)
                    λr2_sol = 1 / (1 / λb_sol - 1 / λr1_fix)
                    push!(pms, CollinearPhaseMatch(cr, (λr1_fix, λr2_sol, λb_sol), temp, hi_or_lo_rrb, θ_ϕ...))
                else
                    fun = λb -> delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r2=λr2_fix, lambda_b=λb, temp=temp)
                    λb_sol = find_zero(fun, (λb_range[i], λb_range[i+1]), Bisection(), atol=tol)
                    λr1_sol = 1 / (1 / λb_sol - 1 / λr2_fix)
                    push!(pms, CollinearPhaseMatch(cr, (λr1_sol, λr2_fix, λb_sol), temp, hi_or_lo_rrb, θ_ϕ...))
                end
            end
        end
    end

    return pms
end

function _pms_vs_lambda_r1(θ_ϕs, hi_or_lo_rrb, cr, λb, temp, λmin, λmax, ngrid, tol)
    pms = CollinearPhaseMatch[]

    for θ_ϕ in θ_ϕs
        if λmax < 1 / (1 / λb - 1 / λmax)
            # if min(λmax, 1 / (1 / λb - 1 / λmax)) < max(λmin, 2 * λb)
            λr1_range = []
        else
            λr1_range = range(
                # max(λmin, 2 * λb),
                # min(λmax, 1 / (1 / λb - 1 / λmax)),
                1 / (1 / λb - 1 / λmax),
                λmax,
                length=ngrid
            )
            Δk = [delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r1=λr1, lambda_b=λb, temp=temp) for λr1 in λr1_range]
        end

        for i in 1:length(λr1_range)-1
            if ustrip(Δk[i] * Δk[i+1]) < 0
                fun = λr1 -> delta_k(θ_ϕ..., hi_or_lo_rrb, cr; lambda_r1=λr1, lambda_b=λb, temp=temp)
                λr1_sol = find_zero(fun, (λr1_range[i], λr1_range[i+1]), Bisection(), atol=tol)
                λr2_sol = 1 / (1 / λb - 1 / λr1_sol)

                all_λ_sol = λr1_sol < λr2_sol ? (λr2_sol, λr1_sol, λb) : (λr1_sol, λr2_sol, λb) # Sort red lambdas, so that λ_r1 is the larger one
                push!(pms, CollinearPhaseMatch(cr, all_λ_sol, temp, hi_or_lo_rrb, θ_ϕ...))
            end
        end
    end

    return pms
end
