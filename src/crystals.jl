export NonlinearCrystal, UnidirectionalCrystal, BidirectionalCrystal, default_lambda, default_temp, is_lambda_valid, valid_lambda_range, optical_axis_angle, RefractionDataHiLo, RefractionData, plot_birefringent_refraction, default_temp

abstract type NonlinearCrystal end

Base.broadcastable(cr::NonlinearCrystal) = Ref(cr)

"""
    default_lambda(cr::NonlinearCrystal)

Returns a default wavelength for use with a crystal `cr`.
If a wavelength range is defined for the X-axis principal index, the midpoint of that range is returned.
Otherwise, the default value is `633 nm`.
"""
function default_lambda(cr::NonlinearCrystal)
    return isnothing(cr.n_X_principal.lambda_range) ? 633u"nm" : sum(cr.n_X_principal.lambda_range) / 2
end

"""
    default_temp(cr::NonlinearCrystal)

Returns the reference temperature from the X-axis principal refractive index of the crystal `cr`.
This is typically used when no explicit temperature is provided.
"""
function default_temp(cr::NonlinearCrystal)
    cr.n_X_principal.temp_ref
end

"""
    is_lambda_valid(lambda::Length, cr::NonlinearCrystal; warn_tol::Length=1u"nm")

Checks whether the given wavelength lies within the valid wavelength range of all three principal refractive index models of the crystal `cr`.
A tolerance `warn_tol` allows small deviations near the boundaries.
Returns `true` if the wavelength is valid for all three axes, otherwise `false`.
"""
function is_lambda_valid(lambda::Length, cr::NonlinearCrystal; warn_tol::Length=1u"nm")
    return is_lambda_valid(lambda, cr.n_X_principal; warn_tol) && is_lambda_valid(lambda, cr.n_Y_principal; warn_tol) && is_lambda_valid(lambda, cr.n_Z_principal; warn_tol)
end

"""
    valid_lambda_range(cr::NonlinearCrystal)

Returns the range of wavelengths that is valid across all three principal axes of the crystal `cr`.
The result is the intersection of the individual wavelength ranges for `n_X`, `n_Y`, and `n_Z`.
"""
function valid_lambda_range(cr::NonlinearCrystal)
    lambda_ranges = [cr.n_X_principal.lambda_range, cr.n_Y_principal.lambda_range, cr.n_Z_principal.lambda_range]
    min_valid_lambda = maximum([l[1] for l in lambda_ranges])
    max_valid_lambda = minimum([l[2] for l in lambda_ranges])
    return (min_valid_lambda, max_valid_lambda)
end

## Unidirectional crystal

"""
    UnidirectionalCrystal(metadata::Dict, n_o_principal::RefractiveIndex,
                          n_e_principal::RefractiveIndex,
                          d_XYZ_ref_full::AbstractArray{<:Number,3};
                          miller_delta=nothing)

Constructs a uniaxial crystal with ordinary (`n_o_principal`) and extraordinary (`n_e_principal`) refractive index models, and a full 3×3×3 nonlinear tensor `d_XYZ_ref_full`.
The `metadata` must include a recognized `:point_group` keyword. If `miller_delta` is provided, it is used for Miller scaling of the effective nonlinear coefficient `d_eff`.
The Z axis is treated as the optical axis.
"""
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
"""
    BidirectionalCrystal(metadata::Dict,
                         n_X_principal::RefractiveIndex,
                         n_Y_principal::RefractiveIndex,
                         n_Z_principal::RefractiveIndex,
                         d_XYZ_ref_full::AbstractArray{<:Number,3};
                         miller_delta=nothing,
                         warn_n_Z_smaller_n_X=true)

Constructs a biaxial crystal using separate refractive index models for the principal X, Y, and Z axes, along with a full nonlinear tensor `d_XYZ_ref_full`.
The `metadata` must include a valid `:point_group`. If `warn_n_Z_smaller_n_X` is true, a warning is issued if the crystal axes are not sorted as `n_Z ≥ n_Y ≥ n_X`, which is the expected convention in this package.
If `miller_delta` is provided, it is used for Miller scaling of the effective nonlinear coefficient `d_eff`.
"""
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

"""
    assign_o_or_e(principal_plane::Symbol, E_dir::AbstractVector{<:Number}; angle_tol_ud::Angle=0.2u"°")

Classifies a polarization vector `E_dir` as ordinary (`:o`) or extraordinary (`:e`) with respect to the given `principal_plane`.
In uniaxial crystals (`:UD`), the classification depends on the angle between `E_dir` and the optical axis. In biaxial crystals, the function assumes `E_dir` lies in or perpendicular to the specified principal plane (`:XY`, `:XZ`, or `:YZ`), and emits an assertion otherwise.
Returns the symbol `:o` or `:e` depending on the geometric configuration.
"""
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

"""
    o_e_to_hi_lo(o_or_e::Symbol, cr::UnidirectionalCrystal, lambda::Length=default_lambda(cr); temp::Temperature=default_temp(cr))

Maps a polarization label `:o` or `:e` to a high (`:hi`) or low (`:lo`) refractive index, based on the comparison between ordinary and extraordinary indices at the given wavelength and temperature.
Returns `:hi` if the refractive index for the given polarization is larger than the other, and `:lo` otherwise.
"""
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

"""
    o_e_to_hi_lo(o_or_e_rrb::NTuple{3,Symbol}, cr::UnidirectionalCrystal, lambda_rrb::NTuple{3,Length}; temp::Temperature=default_temp(cr))

Applies `o_e_to_hi_lo` elementwise to a 3-tuple of polarizations `o_e_to_hi_lo` and a 3-tuple of wavelengths `lambda_rrb`, 
corresponding to the first (typically the longest wavelength) and second red wavelength and the blue wavelength.
Returns a 3-tuple of `:hi`/`:lo` symbols.
"""
function o_e_to_hi_lo(
    o_or_e_rrb::NTuple{3,Symbol},
    cr::UnidirectionalCrystal,
    lambda_rrb::NTuple{3,Length};
    temp::Temperature=default_temp(cr),
)
    return Tuple(o_e_to_hi_lo(o_or_e_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(o_or_e_rrb))
end

"""
    hi_lo_to_o_e(hi_or_lo::Symbol, cr::UnidirectionalCrystal, lambda::Length=default_lambda(cr); temp::Temperature=default_temp(cr))

Maps a `:hi` or `:lo` refractive index label to the corresponding polarization `:o` or `:e`, based on the actual indices at the given wavelength and temperature.
This is the inverse of `o_e_to_hi_lo` and assumes uniaxial behavior.
"""
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

"""
    hi_lo_to_o_e(hi_or_lo_rrb::NTuple{3,Symbol}, cr::UnidirectionalCrystal, lambda_rrb::NTuple{3,Length}; temp::Temperature=default_temp(cr))

Applies `hi_lo_to_o_e` elementwise to a 3-tuple of `:hi`/`:lo` labels and corresponding wavelengths, for use in triple-wave processes.
Returns a tuple of `:o`/`:e` polarization labels.
"""
function hi_lo_to_o_e(
    hi_or_lo_rrb::NTuple{3,Symbol},
    cr::UnidirectionalCrystal,
    lambda_rrb::NTuple{3,Length};
    temp::Temperature=default_temp(cr),
)
    return Tuple(hi_lo_to_o_e(hi_or_lo_rrb[i], cr, lambda_rrb[i]; temp) for i in eachindex(hi_or_lo_rrb))
end

"""
    polarization_rrb_to_hi_lo(pol_rrb::NTuple{3,Symbol}, cr::NonlinearCrystal, lambda_rrb::NTuple{3,Length}; temp::Temperature=default_temp(cr))

Converts a 3-tuple of polarization symbols (`:o`/`:e` or `:hi`/`:lo`) to `:hi`/`:lo`, depending on the type of crystal and the refractive indices at the given wavelengths.
Only uniaxial crystals support conversion from `:o`/`:e`; other crystal types require the input to already use `:hi`/`:lo`.
Returns a 3-tuple of `:hi`/`:lo` symbols or emits an error if the input is invalid.
"""
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

"""
    optical_axis_angle(cr::BidirectionalCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))

Computes the angle between the optical axis and the principal Z-axis for a biaxial crystal `cr`, based on the principal refractive indices at the specified wavelength and temperature.
Returns an angle between 0 and 90°.
"""
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
    # Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#dmitriev2013handbook
    if n_X < n_Z
        Vz = asin((n_Z * sqrt(n_Y^2 - n_X^2)) / (n_Y * sqrt(n_Z^2 - n_X^2)))
    else
        Vz = acos((n_X * sqrt(n_Y^2 - n_Z^2)) / (n_Y * sqrt(n_X^2 - n_Z^2)))
    end

    return Vz |> u"rad"
end

"""
    optical_axis_angle(cr::UnidirectionalCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))

Returns the optical axis angle for a uniaxial crystal, which is always zero by definition since the optical axis aligns with the principal Z-axis.
"""
function optical_axis_angle(
    cr::UnidirectionalCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr),
)
    return 0.0u"rad"
end

"""
    RefractionTypeHiLo

Describes the polarization type of the two birefringent solutions (`:hi` and `:lo`) in a given principal plane.
Each polarization is labeled as ordinary (`:o`) or extraordinary (`:e`) based on field orientation.
"""
struct RefractionTypeHiLo
    principal_plane::Union{Symbol}
    o_or_e_hi_lo::NTuple{2,Symbol}
end

"""
    RefractionTypeHiLo(principal_plane::Symbol, E_dir_hi_lo::NTuple{2,<:AbstractVector{<:Number}})

Classifies a pair of electric field directions as ordinary or extraordinary relative to the specified principal plane. 
Used to track the polarization type of high and low index branches in birefringent crystals.
"""
function RefractionTypeHiLo(
    principal_plane::Symbol,
    E_dir_hi_lo::NTuple{2,<:AbstractVector{<:Number}};
)
    o_or_e_hi_lo = Tuple(assign_o_or_e(principal_plane, E_dir) for E_dir in E_dir_hi_lo)

    return RefractionTypeHiLo(principal_plane, o_or_e_hi_lo)
end

"""
    RefractionDataHiLo

Contains all relevant optical data for the two birefringent branches (`:hi` and `:lo`) in a nonlinear crystal.
Includes refractive indices, group indices, energy flow vectors, walk-off angles, and dispersion derivatives (β₀ through β₃).
Computed for a specific propagation direction, wavelength, and temperature.
`RefractionDataHiLo` can be split into [`RefractionData`](@ref) instances to represent only one of both polarization branches.
"""
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

    @printf(io, "%-25s %-25s %-25s\n", "Refractive index:", auto_fmt(rd_hi.n; digits), auto_fmt(rd_lo.n; digits))
    @printf(io, "%-25s %-25s %-25s\n", "Group index:", auto_fmt(rd_hi.group_index), auto_fmt(rd_lo.group_index))
    @printf(io, "%-25s %-25s %-25s\n", "Walkoff angle (mrad):",
        auto_fmt(ustrip(u"mrad", rd_hi.walkoff_angle); digits),
        auto_fmt(ustrip(u"mrad", rd_lo.walkoff_angle); digits),
    )
    @printf(io, "%-25s %-25s %-25s\n", "S direction:", vec_str(rd_hi.S_dir), vec_str(rd_lo.S_dir))
    @printf(io, "%-24s %-25s %-25s\n", "E direction:", "±" * vec_str(rd_hi.E_dir), "±" * vec_str(rd_lo.E_dir))
    @printf(io, "%-24s %-25s %-25s\n", "D direction:", "±" * vec_str(rd_hi.D_dir), "±" * vec_str(rd_lo.D_dir))
    @printf(io, "%-25s %-25s %-25s\n", "β₂ (fs²/mm):",
        auto_fmt(ustrip(u"fs^2/mm", rd_hi.beta2); digits),
        auto_fmt(ustrip(u"fs^2/mm", rd_lo.beta2); digits),
    )
    @printf(io, "%-25s %-25s %-25s\n", "β₃ (fs³/mm):",
        auto_fmt(ustrip(u"fs^3/mm", rd_hi.beta3); digits),
        auto_fmt(ustrip(u"fs^3/mm", rd_lo.beta3); digits),
    )
end

"""
    RefractionDataHiLo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr);
                       temp::Temperature=default_temp(cr), angle_tol::Angle=0.1u"°")

Constructs a high/low index birefringent refraction object `RefractionDataHiLo` for a given propagation direction defined by spherical angles `theta`, `phi`. 
Internally calculates refractive indices, energy and phase directions, walk-off angles, group indices, and frequency derivatives 
for both polarization branches. Polarization classification is automatic based on principal planes or optical axis orientation.
"""
function RefractionDataHiLo(
    theta::Angle,
    phi::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
    angle_tol::Angle=0.1u"°",
)
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(theta, phi, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)

    E_dir_hi_lo, S_dir_hi_lo = calc_E_dir_S_dir(k_dir, ε_tensor, D_dir_hi_lo)

    if isa(cr, UnidirectionalCrystal)
        # o and e are valid classifiers for UnidirectionalCrystal even if propagating outside of the pricipal planes
        refr_type_hi_lo = (RefractionTypeHiLo(:UD, E_dir_hi_lo), nothing)
    else
        principal_planes = find_principal_planes(theta, phi; angle_tol)
        refr_type_hi_lo = Tuple((isnothing(p) ? nothing : RefractionTypeHiLo(p, E_dir_hi_lo) for p in principal_planes))
    end

    walkoff_angle_hi_lo = (s -> (acos(clamp(dot(s, k_dir), -1, 1)))).(S_dir_hi_lo) .|> u"rad"

    # Calculate derivative based data
    group_index_hi_lo = calc_group_index_hi_lo(theta, phi, cr, lambda, temp)
    β0_hi_lo = calc_β0_hi_lo(theta, phi, cr, lambda, temp)
    β1_hi_lo = calc_β1_hi_lo(theta, phi, cr, lambda, temp)
    β2_hi_lo = calc_β2_hi_lo(theta, phi, cr, lambda, temp)
    β3_hi_lo = calc_β3_hi_lo(theta, phi, cr, lambda, temp)

    rd = RefractionDataHiLo{typeof(cr)}(
        theta,
        phi,
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

"""
    calc_n_hi_lo(θ::Angle, ϕ::Angle, cr::NonlinearCrystal, lambda::Length; temp::Temperature)

Computes the refractive indices of the two birefringent eigenmodes for a given direction in spherical angles. 
This is a minimal calculation used for derivative tracing (e.g., with ForwardDiff) and does not return field vectors.
"""
function calc_n_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length; temp::Temperature)
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(theta, phi, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    return n_hi_lo
end

"""
    calc_k_dir_ε_tensor_n_XYZ(θ::Angle, ϕ::Angle, cr::NonlinearCrystal,
                              lambda::Length=default_lambda(cr),
                              temp::Temperature=default_temp(cr))

Computes the wavevector direction from spherical angles and constructs the dielectric tensor ε for the given crystal `cr` and conditions. 
Also returns the principal refractive indices along the X, Y, and Z axes for the given wavelength `lambda` and temperature `temp`.
"""
function calc_k_dir_ε_tensor_n_XYZ(
    theta::Angle,
    phi::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr),
)
    n_X = refractive_index(cr.n_X_principal, lambda, temp)
    n_Y = refractive_index(cr.n_Y_principal, lambda, temp)
    n_Z = refractive_index(cr.n_Z_principal, lambda, temp)

    k_dir = angles_to_vector(theta, phi)
    ε_tensor = Diagonal(@SVector [1 / n_X^2, 1 / n_Y^2, 1 / n_Z^2])

    return k_dir, ε_tensor, (n_X, n_Y, n_Z)
end

"""
    calc_n_hi_lo_D_dir_hi_lo(k_dir::AbstractVector{<:Real}, ε_tensor::AbstractMatrix{<:Real})

Given a propagation direction `k_dir` and dielectric tensor `ε_tensor`, computes the two effective refractive indices and corresponding 
displacement vectors (`D`) for the birefringent modes. Works by diagonalizing the projected dielectric tensor in the plane orthogonal to `k_dir`.
"""
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

"""
    calc_E_dir_S_dir(k_dir::AbstractVector{<:Real},
                     ε_tensor::AbstractMatrix{<:Real},
                     D_dir_hi_lo::Tuple)

From the `D` vectors and dielectric tensor, computes the electric field directions `E = ε⁻¹·D` and corresponding 
Poynting vectors `S = E × H`, where `H = k × E`. The resulting vectors are normalized and aligned with `k_dir`.
"""
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
function calc_β0_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    k_dir, ε_tensor, n_XYZ = calc_k_dir_ε_tensor_n_XYZ(theta, phi, cr, lambda, temp)
    n_hi_lo, D_dir_hi_lo = calc_n_hi_lo_D_dir_hi_lo(k_dir, ε_tensor)
    return Tuple([(2π / lambda * n) |> u"m^-1" for n in n_hi_lo])
end
calc_β0_hi_lo(theta::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β0_hi_lo(theta, 0.0u"rad", cr, args...; kwargs...)

# β1 = ∂k/∂ω
function calc_β1_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β0_hi_lo(theta, phi, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s/m")
end
calc_β1_hi_lo(theta::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β1_hi_lo(theta, 0.0u"rad", cr, args...; kwargs...)

# β2 = ∂²k/∂ω²
function calc_β2_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β1_hi_lo(theta, phi, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^2/m")
end
calc_β2_hi_lo(theta::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β2_hi_lo(theta, 0.0u"rad", cr, args...; kwargs...)

# β3 = ∂³k/∂ω³
function calc_β3_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    ω_in = 2π * c_0 / lambda

    fun = ω -> [ustrip(d) for d in calc_β2_hi_lo(theta, phi, cr, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)]
    return Tuple(
        ForwardDiff.derivative(
            fun,
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^3/m")
end
calc_β3_hi_lo(theta::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_β3_hi_lo(theta, 0.0u"rad", cr, args...; kwargs...)

function calc_group_index_hi_lo(theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr), temp::Temperature=default_temp(cr))
    return Tuple([β1 * c_0 for β1 in calc_β1_hi_lo(theta, phi, cr, lambda, temp)])
end
calc_group_index_hi_lo(theta::Angle, cr::UnidirectionalCrystal, args...; kwargs...) = calc_group_index_hi_lo(theta, 0.0u"rad", cr, args...; kwargs...)

"""
    plot_birefringent_refraction(theta, phi, cr::NonlinearCrystal,
                                  lambda::Length=default_lambda(cr),
                                  temp::Temperature=default_temp(cr))

Creates a 3D GLMakie visualization of birefringent refraction for the given direction and crystal. The plot shows the index ellipsoid, 
wavevector direction `k`, and the Poynting vector (`S`), electric fields (`E`), and displacement (`D`) vectors for both high and low index solutions. 
Useful for verifying polarization, walk-off behavior, and eigenvector orientation.
"""
function plot_birefringent_refraction(
    theta,
    phi,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr),
    temp::Temperature=default_temp(cr);
)
    k_dir, ε_tensor, (n_X, n_Y, n_Z) = calc_k_dir_ε_tensor_n_XYZ(theta, phi, cr, lambda, temp)
    rd = RefractionDataHiLo(theta, phi, cr, lambda; temp)

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

"""
    RefractionType

Describes a single polarization state (`:o` or `:e`) in a specified principal plane.
"""
struct RefractionType
    principal_plane::Union{Symbol}
    o_or_e::Symbol
end

"""
    RefractionType(hi_or_lo::Symbol, rt::RefractionTypeHiLo)

Extracts a single-polarization `RefractionType` from a two-branch `RefractionTypeHiLo`, corresponding to either the high or low index solution.
"""
function RefractionType(hi_or_lo::Symbol, rt::RefractionTypeHiLo)
    return RefractionType(rt.principal_plane, rt.o_or_e_hi_lo[hi_or_lo === :hi ? 1 : 2])
end

"""
    RefractionData

A single-branch slice of [`RefractionDataHiLo`](@ref) corresponding to either the `:hi` or `:lo` polarization branch.
Stores all relevant scalar and vector optical quantities for one mode: refractive index, group index, Poynting vector, polarization directions, walk-off, and β dispersion parameters.
"""
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

    @printf(io, "%-25s %-25s\n", "Refractive index:", auto_fmt(rd.n; digits))
    @printf(io, "%-25s %-25s\n", "Group index:", auto_fmt(rd.group_index))
    @printf(io, "%-25s %-25s\n", "Walkoff angle (mrad):",
        auto_fmt(ustrip(u"mrad", rd.walkoff_angle); digits),
    )
    @printf(io, "%-25s %-25s\n", "S direction:", vec_str(rd.S_dir))
    @printf(io, "%-24s %-25s\n", "E direction:", "±" * vec_str(rd.E_dir))
    @printf(io, "%-24s %-25s\n", "D direction:", "±" * vec_str(rd.D_dir))
    @printf(io, "%-25s %-25s\n", "β₂ (fs²/mm):",
        auto_fmt(ustrip(u"fs^2/mm", rd.beta2); digits),
    )
    @printf(io, "%-25s %-25s\n", "β₃ (fs³/mm):",
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

"""
    RefractionData(hi_or_lo::Symbol, theta::Angle, phi::Angle, cr::NonlinearCrystal, lambda::Length=default_lambda(cr);
                   temp::Temperature=default_temp(cr))

Convenience constructor that combines computation and extraction: builds a `RefractionDataHiLo` and returns the 
single-polarization `RefractionData` corresponding to `:hi` or `:lo`.
"""
function RefractionData(
    hi_or_lo::Symbol,
    theta::Angle,
    phi::Angle,
    cr::NonlinearCrystal,
    lambda::Length=default_lambda(cr);
    temp::Temperature=default_temp(cr),
)
    return RefractionData(hi_or_lo, RefractionDataHiLo(theta, phi, cr, lambda; temp))
end

"""
    plot_refractiveindex(cr::BidirectionalCrystal; n_sample_pts=500,
                         temp=[default_temp(cr)], lambda_min=nothing, lambda_max=nothing)

Plots the refractive indices `n_X`, `n_Y`, and `n_Z` as functions of wavelength for a biaxial crystal using GLMakie.
Each principal axis is drawn with a different color and can optionally show temperature dependence.
Returns a figure with legend and interactive inspection enabled.
"""
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

"""
    plot_refractiveindex(cr::UnidirectionalCrystal; n_sample_pts=500,
                         temp=[default_temp(cr)], lambda_min=nothing, lambda_max=nothing)

Plots the ordinary (`n_X = n_Y`) and extraordinary (`n_Z`) refractive indices of a uniaxial crystal over wavelength using GLMakie.
The ordinary index is shown once for both `X` and `Y`, and different colormaps distinguish the two polarizations.
Returns a Makie figure with labeled axes, legend, and interactive inspection.
"""
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