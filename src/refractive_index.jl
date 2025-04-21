export SellmeierFunction, refractive_index, group_index, phase_velocity, group_velocity, β0, β1, β2, β3, dn_dtemp, plot_refractiveindex, plot_refractiveindex!

# Type piracy: Solve ambiguity during automatic differentiation with units # TODO: Is this always correct? Is there a cleaner way?
Base.convert(::Type{ForwardDiff.Dual{T,V,N}}, x::Quantity) where {N,V,T} = uconvert(Unitful.NoUnits, x)

abstract type RefractiveIndex end

Base.broadcastable(ri::RefractiveIndex) = Ref(ri)

(ri::RefractiveIndex)(lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri)) = refractive_index(ri, lambda, temp) # Shorthand for refractive_index

"""
    default_lambda(ri::RefractiveIndex)

Returns a default wavelength for use with the refractive index model `ri`.
If `ri.lambda_range` is defined, it returns the midpoint of the range; otherwise, it defaults to 633 nm.
"""
function default_lambda(ri::RefractiveIndex)
    return isnothing(ri.lambda_range) ? 633u"nm" : sum(ri.lambda_range) / 2
end

"""
    default_temp(ri::RefractiveIndex)

Returns the reference temperature `temp_ref` defined in the refractive index model `ri`.
This is typically used as the default when no explicit temperature is provided.
"""
function default_temp(ri::RefractiveIndex)
    ri.temp_ref
end

# Compute first, second, and third derivatives using ForwardDiff.jl
function dn_dλ(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return ForwardDiff.derivative(lambda -> refractive_index(ri, lambda * 1u"m", temp), ustrip(lambda |> u"m")) * 1u"m^-1"
end

function d2n_dλ2(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return ForwardDiff.derivative(lambda -> ustrip(dn_dλ(ri, lambda * 1u"m", temp)), ustrip(lambda |> u"m")) * 1u"m^-2"
end

function d3n_dλ3(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return ForwardDiff.derivative(lambda -> ustrip(d2n_dλ2(ri, lambda * 1u"m", temp)), ustrip(lambda |> u"m")) * 1u"m^-3"
end

"""
    β0(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Returns the wavevector `β₀ = 2π/λ · n(λ)` for a given wavelength and temperature, representing the phase accumulation per meter.
The result has units of `1/m`.
"""
# β0 = k
function β0(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return (2π / lambda * refractive_index(ri, lambda, temp)) |> u"m^-1"
end

"""
    β1(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Calculates the first derivative of the wavevector `β` with respect to angular frequency `ω`.
This value determines the group delay and contributes to dispersion calculations.
"""
# β1 = ∂k/∂ω
function β1(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β0(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s/m")
end

"""
    β2(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Computes the second-order dispersion parameter by differentiating `β₁` with respect to angular frequency.
This characterizes group velocity dispersion (GVD), e.g. important for ultrafast pulse propagation.
"""
# β2 = ∂²k/∂ω²
function β2(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β1(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^2/m")
end

"""
    β3(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Calculates the third-order derivative of the wavevector with respect to angular frequency.
This quantifies third-order dispersion (TOD), e.g. relevant for modeling higher-order effects in broadband pulse propagation.
"""
# β3 = ∂³k/∂ω³
function β3(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β2(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^3/m")
end

"""
    group_index(ri, lambda=default_lambda(ri), temp=default_temp(ri))

Computes the group index `n_g`, which determines the group velocity in a dispersive medium.
This is calculated as `β₁ · c`, where `β₁` is the first derivative of the propagation constant with respect to angular frequency.
The result is unitless.
"""
function group_index(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return β1(ri, lambda, temp) * c_0
end

"""
    phase_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Returns the phase velocity `v_p = c / n(λ)` of a monochromatic wave in the material described by `ri` at the given conditions.
The result has physical units of velocity.
"""
function phase_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return c_0 / refractive_index(ri, lambda, temp) |> u"m/s"
end

"""
    group_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Computes the group velocity of a wavepacket centered at `lambda`, using the first derivative of the wavevector with respect to frequency.
This quantity reflects how the envelope of a pulse propagates through the medium.
"""
function group_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return 1 / β1(ri, lambda, temp) |> u"m/s"
end

"""
    dn_dtemp(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))

Computes the temperature derivative of the `RefractiveIndex` at the specified wavelength `lambda` and temperature `temp`.
This is used to analyze thermal drift, tuning, or thermal lensing effects in nonlinear and resonant systems.
"""
function dn_dtemp(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return ForwardDiff.derivative(
        t -> ustrip(refractive_index(ri, lambda, t * 1u"K")),
        ustrip(u"K", temp)
    ) * 1u"K^-1"
end

"""
    SellmeierFunction(n_fun, lambda_range=nothing; temp_ref=293.15u"K", temp_range=nothing)

Wraps a wavelength- and temperature-dependent refractive index function in a `RefractiveIndex` model. 
This is typically used to represent Sellmeier equations or any function that evaluates `n(λ, T)`, possibly dependent on a reference temperature `temp_ref`.
If `lambda_range` and `temp_range` are provided, they define the validity range and are checked optionally during evaluation.
"""
struct SellmeierFunction{FL} <: RefractiveIndex
    n_fun::FL
    lambda_range::NTuple{2,typeof(1.0u"m")}
    temp_ref::typeof(1.0u"K")
    temp_range::Union{Nothing,NTuple{2,typeof(1.0u"K")}}

    function SellmeierFunction(
        n_fun::Function,
        lambda_range::Union{Nothing,NTuple{2,Length}}=nothing;
        temp_ref=293.15u"K",
        temp_range::Union{Nothing,NTuple{2,Temperature}}=nothing,
    )
        temp_ref = temp_ref |> u"K"
        return new{typeof(n_fun)}(
            n_fun,
            lambda_range,
            temp_ref,
            temp_range
        )
    end
end

"""
    is_lambda_valid(lambda::Length, sri::SellmeierFunction; warn_tol::Length=1u"nm")

Checks whether the given wavelength is within the valid `lambda_range` of a `SellmeierFunction`.
A tolerance `warn_tol` can be specified to allow for small numerical deviations near the boundaries.
Returns `true` if the wavelength lies within the specified valid range, `false` otherwise.
"""
function is_lambda_valid(lambda::Length, sri::SellmeierFunction; warn_tol::Length=1u"nm")
    return lambda ≥ (sri.lambda_range[1] - warn_tol) && lambda ≤ (sri.lambda_range[2] + warn_tol)
end

"""
    refractive_index(sri::SellmeierFunction, lambda::Length, temp::Temperature=sri.temp_ref;
                     check_lambda_range::Symbol=:warn, check_temp_range::Symbol=:warn,
                     warn_tol::Length=1u"nm")

Evaluates the refractive index for a given wavelength and temperature using the stored function inside a `SellmeierFunction`.
If bounds are set and `check_lambda_range` or `check_temp_range` is not `:none`, warnings or errors are issued for out-of-range inputs.
Returns a unitless number corresponding to the refractive index.
"""
function refractive_index(
    sri::SellmeierFunction,
    lambda::Length,
    temp::Temperature=sri.temp_ref;
    check_lambda_range::Symbol=:warn,
    check_temp_range::Symbol=:warn,
    warn_tol::Length=1u"nm",
)
    if check_lambda_range != :none
        if !is_lambda_valid(lambda, sri; warn_tol)
            str = "Wavelength λ=$(lambda) out of valid range for Sellmeier equation [$(sri.lambda_range[1]), $(sri.lambda_range[2])]"
            check_lambda_range == :warn ? (@warn str) : (@error str)
        end
    end

    if !isnothing(sri.temp_range) && check_temp_range != :none
        if lambda < sri.lambda_range[1] || lambda > sri.lambda_range[2]
            str = "Temperature T=$(temp) out of valid range for Sellmeier equation [$(sri.temp_range[1]), $(sri.temp_range[2])]"
            check_temp_range == :warn ? (@warn str) : (@error str)
        end
    end

    return uconvert(Unitful.NoUnits, sri.n_fun(lambda, temp))
end

## Plots

"""
    plot_refractiveindex(ri::RefractiveIndex; n_sample_pts=500, temp=nothing, lambda_min=nothing, lambda_max=nothing, label="", colormap=:vik)

Creates and returns a new `GLMakie` plot of the refractive index as a function of wavelength, optionally at multiple temperatures if 
    a vector of temperatures is provided. If `temp` is omitted, the reference temperature of `ri` is used. 
"""
function plot_refractiveindex(
    ri::RefractiveIndex;
    n_sample_pts=500,
    temp=nothing,
    lambda_min=nothing,
    lambda_max=nothing,
    label="",
    colormap=:vik,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")
    plot_refractiveindex!(ri; n_sample_pts, temp, lambda_min, lambda_max, label, colormap)
    Legend(f[1, 2], ax)
    return f
end

"""
    plot_refractiveindex!(ri::RefractiveIndex; n_sample_pts=500, temp=nothing, lambda_min=nothing,
                          lambda_max=nothing, label="", colormap=:vik, digits::Integer=3)

Adds a refractive index vs. wavelength curve to the current `GLMakie` axis, optionally for multiple temperatures.
This is useful for overlaying plots from different models or conditions. The interactive inspector shows wavelength, index, and temperature.
"""
function plot_refractiveindex!(
    ri::RefractiveIndex;
    n_sample_pts=500,
    temp=nothing,
    lambda_min=nothing,
    lambda_max=nothing,
    label="",
    colormap=:vik,
    digits::Integer=3,
)
    temp = isnothing(temp) ? ri.temp_ref : temp
    lambda_min = isnothing(lambda_min) ? ri.lambda_range[1] : lambda_min
    lambda_max = isnothing(lambda_max) ? ri.lambda_range[2] : lambda_max
    temp = isnothing(temp) ? ri.temp_ref : temp
    λs = LinRange(lambda_min, lambda_max, n_sample_pts)

    isempty(label) || (label = label * ": ")
    for (i, Ti) in enumerate(temp)
        ris = refractive_index.(ri, λs, Ti)
        inspector_label = (plot, index, position) -> label * "\nλ=$(round(u"µm", position[1] * u"µm"; digits))\nn=$(round(position[2]; digits))\nT = $(Ti)"
        lines!(λs, ris; label=label * "T = $(Ti)", color=i, colormap, colorrange=(1, length(temp) + 2), inspectable=true, inspector_label)
    end
end