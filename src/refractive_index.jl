export SellmeierFunction, refractive_index, group_index, phase_velocity, group_velocity, β0, β1, β2, β3, dn_dtemp, plot_refractiveindex, plot_refractiveindex!

# Type piracy: Solve ambiguity during automatic differentiation with units # TODO: Is this always correct? Is there a cleaner way?
Base.convert(::Type{ForwardDiff.Dual{T,V,N}}, x::Quantity) where {N,V,T} = uconvert(Unitful.NoUnits, x)

abstract type RefractiveIndex end

Base.broadcastable(ri::RefractiveIndex) = Ref(ri)

(ri::RefractiveIndex)(lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri)) = refractive_index(ri, lambda, temp) # Shorthand for refractive_index

function default_lambda(ri::RefractiveIndex)
    return isnothing(ri.lambda_range) ? 633u"nm" : sum(ri.lambda_range) / 2
end

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

# β0 = k
function β0(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return (2π / lambda * refractive_index(ri, lambda, temp)) |> u"m^-1"
end

# β1 = ∂k/∂ω
function β1(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β0(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s/m")
end

# β2 = ∂²k/∂ω²
function β2(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β1(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^2/m")
end

# β3 = ∂³k/∂ω³
function β3(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    ω_in = 2π * c_0 / lambda
    return (
        ForwardDiff.derivative(
            ω -> ustrip(β2(ri, uconvert(u"m", 2π * c_0 / ω * 1u"s"), temp)),
            ustrip(ω_in |> u"s^-1")
        ) * 1u"s^3/m")
end


function group_index(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return β1(ri, lambda, temp) * c_0
end

function phase_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return c_0 / refractive_index(ri, lambda, temp) |> u"m/s"
end

function group_velocity(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return 1 / β1(ri, lambda, temp) |> u"m/s"
end

function dn_dtemp(ri::RefractiveIndex, lambda::Length=default_lambda(ri), temp::Temperature=default_temp(ri))
    return ForwardDiff.derivative(
        t -> ustrip(refractive_index(ri, lambda, t * 1u"K")),
        ustrip(u"K", temp)
    ) * 1u"K^-1"
end


struct SellmeierFunction{FL,LR,TF,TR} <: RefractiveIndex
    n_fun::FL
    lambda_range::LR
    temp_ref::TF
    temp_range::TR

    function SellmeierFunction(
        n_fun::Function,
        lambda_range::Union{Nothing,Tuple{Unitful.Length,Unitful.Length}}=nothing;
        temp_ref=293.15u"K",
        temp_range::Union{Nothing,Tuple{Unitful.Temperature,Unitful.Temperature}}=nothing,
    )
        temp_ref = temp_ref |> u"K"
        return new{eltype(n_fun),typeof(lambda_range),typeof(temp_ref),typeof(temp_range)}(
            n_fun,
            lambda_range,
            temp_ref,
            temp_range
        )
    end
end

function is_lambda_valid(lambda::Length, sri::SellmeierFunction; warn_tol::Length=1u"nm")
    return lambda ≥ (sri.lambda_range[1] - warn_tol) && lambda ≤ (sri.lambda_range[2] + warn_tol)
end

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