export SellmeierFunction, refractive_index, group_index, phase_velocity, group_velocity, β0, β1, β2, β3, group_velocity_dispersion, third_order_dispersion, plot_refractiveindex, plot_refractiveindex!

# Solve ambiguity during automatic differentiation with units # TODO: Is this always correct?
Base.convert(::Type{ForwardDiff.Dual{T,V,N}}, x::Quantity) where {N,V,T} = uconvert(Unitful.NoUnits, x)

abstract type RefractiveIndex end

Base.broadcastable(ri::RefractiveIndex) = Ref(ri)

(ri::RefractiveIndex)(λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref) = refractive_index(ri, λ, T) # Shorthand for refractive_index

# Compute first, second, and third derivatives using ForwardDiff.jl
function dn_dλ(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return ForwardDiff.derivative(λ -> refractive_index(ri, λ * 1u"m", T), ustrip(λ |> u"m")) * 1u"m^-1"
end

function d2n_dλ2(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return ForwardDiff.derivative(λ -> ustrip(dn_dλ(ri, λ * 1u"m", T)), ustrip(λ |> u"m")) * 1u"m^-2"
end

function d3n_dλ3(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return ForwardDiff.derivative(λ -> ustrip(d2n_dλ2(ri, λ * 1u"m", T)), ustrip(λ |> u"m")) * 1u"m^-3"
end

# β0 = k
function β0(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return (2π / λ * refractive_index(ri, λ, T)) |> u"m^-1"
end

# β1 = ∂k/∂ω
function β1(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    ω_in = 2π * c_0 / λ
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β0(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s/m")
end

# β2 = ∂²k/∂ω²
function β2(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    ω_in = 2π * c_0 / λ
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β1(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^2/m")
end

# β3 = ∂³k/∂ω³
function β3(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    ω_in = 2π * c_0 / λ
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β2(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), T)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^3/m")
end


function group_index(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return β1(ri, λ, T) * c_0
end

function phase_velocity(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return c_0 / refractive_index(ri, λ, T) |> u"m/s"
end

function group_velocity(ri::RefractiveIndex, λ::Unitful.Length, T::Unitful.Temperature=ri.T_ref)
    return 1 / β1(ri, λ, T) |> u"m/s"
end


# Sellmeier equation: n = f1(λ,T), dn_dT = f2(λ,T)
# Functions f1 and f2 must be selected to match the final evaluation of the refractive index according to: n(λ,T) = f1(λ,T) + f2(λ,T) * (T - T_ref) 
struct SellmeierFunction{FL,LR,TF,FT,TR} <: RefractiveIndex
    n_fun::FL
    lambda_range::LR
    T_ref::TF
    dn_dT_fun::FT
    T_range::TR

    function SellmeierFunction(
        n_fun::Function,
        lambda_range::Union{Nothing,Tuple{Unitful.Length,Unitful.Length}}=nothing;
        T_ref=293.15u"K",
        dn_dT_fun::Function=(λ, T) -> 0.0u"K^-1",
        T_range::Union{Nothing,Tuple{Unitful.Temperature,Unitful.Temperature}}=nothing,
    )
        T_ref = T_ref |> u"K"
        return new{eltype(n_fun),typeof(lambda_range),typeof(T_ref),typeof(dn_dT_fun),typeof(T_range)}(
            n_fun,
            lambda_range,
            T_ref,
            dn_dT_fun,
            T_range
        )
    end
end

function refractive_index(
    sri::SellmeierFunction,
    λ::Unitful.Length,
    T::Unitful.Temperature=sri.T_ref;
    check_lambda_range::Symbol=:warn,
    check_T_range::Symbol=:warn
)
    if check_lambda_range != :none
        if λ < sri.lambda_range[1] || λ > sri.lambda_range[2]
            str = "Wavelength λ=$(λ) out of valid range for Sellmeier equation [$(sri.lambda_range[1]), $(sri.lambda_range[2])]"
            check_lambda_range == :warn ? (@warn str) : (@error str)
        end
    end

    if !isnothing(sri.T_range) && check_T_range != :none
        if λ < sri.lambda_range[1] || λ > sri.lambda_range[2]
            str = "Temperature T=$(T) out of valid range for Sellmeier equation [$(sri.T_range[1]), $(sri.T_range[2])]"
            check_T_range == :warn ? (@warn str) : (@error str)
        end
    end

    return uconvert(Unitful.NoUnits, sri.n_fun(λ, T) + sri.dn_dT_fun(λ, T) * (T - sri.T_ref))
end

## Plots

function plot_refractiveindex(
    ri::RefractiveIndex;
    n_sample_pts=500,
    T=nothing,
    lambda_min=nothing,
    lambda_max=nothing,
    label=nothing,
    colormap=nothing,
)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Wavelength", ylabel="Refractive index")
    plot_refractiveindex!(ri; n_sample_pts, T, lambda_min, lambda_max, label, colormap)
    Legend(f[1, 2], ax)
    return f
end

# TODO: Add option for legend and color and label
function plot_refractiveindex!(
    ri::RefractiveIndex;
    n_sample_pts=500,
    T=nothing,
    lambda_min=nothing,
    lambda_max=nothing,
    label=nothing,
    colormap=nothing,
    digits::Integer=3,
)
    T = isnothing(T) ? ri.T_ref : T
    lambda_min = isnothing(lambda_min) ? ri.lambda_range[1] : lambda_min
    lambda_max = isnothing(lambda_max) ? ri.lambda_range[2] : lambda_max
    T = isnothing(T) ? ri.T_ref : T
    λs = LinRange(lambda_min, lambda_max, n_sample_pts)

    for (i, Ti) in enumerate(T)
        ris = refractive_index.(ri, λs, Ti)
        inspector_label = (plot, index, position) -> label * "\nλ=$(round(u"µm", position[1] * u"µm"; digits))\nn=$(round(position[2]; digits))\nT = $(Ti)"
        lines!(λs, ris; label=label * ": T = $(Ti)", color=i, colormap, colorrange=(1, length(T)+2), inspectable=true, inspector_label)
    end
end