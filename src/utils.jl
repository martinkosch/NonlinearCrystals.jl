export draw_ellipsoid!

function draw_ellipsoid!(
    ax,
    center::AbstractVector=[0.0, 0.0, 0.0],
    principal_axes::AbstractVector=[1.0, 1.0, 1.0];
    # rotation::AbstractVector=[0.0, 0.0, 0.0];
    kwargs...,
)
    s = Sphere(Point3f(center...), 1f0)

    meshplot = mesh!(ax, s, alpha=0.3, transparency=true)
    scale!(meshplot, principal_axes...)

    return ax
end

function vector_to_angles(v::AbstractVector{<:Real})
    x, y, z = v
    r = sqrt(x^2 + y^2 + z^2)
    θ = Float64(acos(z / r)) |> u"rad"   # Polar angle (0 to π)
    ϕ = atan(y, x)  # Azimuthal angle (−π to π)
    ϕ = mod(ϕ, 2π) |> u"rad"    # Normalize φ to [0, 2π)
    return θ, ϕ
end

function angles_to_vector(θ, ϕ)
    x = sin(θ) * cos(ϕ)
    y = sin(θ) * sin(ϕ)
    z = cos(θ)
    return @SVector [x, y, z]
end