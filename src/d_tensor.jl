export calc_d_XYZ_full, rot_mat_crys_to_diel, plot_axes_assignment_crys_to_diel, calc_miller_delta

"""
    expand_voigt_index(i::Integer, l::Integer) -> Vector{NTuple{3,Int}}

Given contracted Voigt indices (i, l), return list of (i, j, k) index tuples.
Some `l` values map to two symmetric (j, k) combinations.
"""
function expand_voigt_index(i::Integer, l::Integer)::Vector{NTuple{3,Int}}
    # Map contracted index l to full (j, k) pairs
    jk_pairs = Dict(
        1 => [(1, 1)],
        2 => [(2, 2)],
        3 => [(3, 3)],
        4 => [(2, 3), (3, 2)],
        5 => [(1, 3), (3, 1)],
        6 => [(1, 2), (2, 1)]
    )

    return [(i, j, k) for (j, k) in jk_pairs[l]]
end

"""
    contract_voigt_index(i::Integer, j::Integer, k::Integer) -> (i::Int, l::Int)

Given full tensor indices (i, j, k), return the contracted Voigt index (i, l).
This assumes symmetry in (j, k): for example, (1,3) and (3,1) both map to l = 5.
"""
function contract_voigt_index(i::Integer, j::Integer, k::Integer)::Tuple{Int,Int}
    # Map sorted (j,k) → l
    jk_to_l = Dict(
        (1, 1) => 1,
        (2, 2) => 2,
        (3, 3) => 3,
        (2, 3) => 4, (3, 2) => 4,
        (1, 3) => 5, (3, 1) => 5,
        (1, 2) => 6, (2, 1) => 6
    )

    l = jk_to_l[(j, k)]
    return (i, l)
end

function tensor_indices(comp::Symbol)
    s = string(comp)
    length(s) == 3 && s[1] == 'd' || error("Invalid tensor component symbol: $comp")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    (i, j)
end

# Expand 3x6 contracted tensor into full 3x3x3
function expand_d_contract(d_contract::AbstractMatrix{<:Number})
    d_full = zeros(eltype(d_contract), 3, 3, 3)
    idx_map = [
        (1, 1), (2, 2), (3, 3),  # J=1,2,3
        (2, 3), (1, 3), (1, 2)   # J=4,5,6 (symmetric indices)
    ]
    for i in 1:3, J in 1:6
        j, k = idx_map[J]
        d_full[i, j, k] += d_contract[i, J]
        if j != k
            d_full[i, k, j] += d_contract[i, J]  # symmetry
        end
    end
    return d_full
end

# Rotate 3x3x3 tensor
function rotate_tensor3(d::AbstractArray{<:Number,3}, R::AbstractMatrix{<:Number})
    @tullio d_rot[i_dash, j_dash, k_dash] := R[i_dash, i] * R[j_dash, j] * R[k_dash, k] * d[i, j, k]
end

# Contract back to 3x6
function contract_d_tensor(d_full::AbstractArray{<:Number,3})
    d_contract = zeros(eltype(d_full), 3, 6)
    idx_map = [
        (1, 1), (2, 2), (3, 3),
        (2, 3), (1, 3), (1, 2)
    ]
    for i in 1:3, J in 1:6
        j, k = idx_map[J]
        if j == k
            d_contract[i, J] = d_full[i, j, k]
        else
            d_contract[i, J] = 0.5 * (d_full[i, j, k] + d_full[i, k, j])
        end
    end
    return d_contract
end

"""
    rot_mat_crys_to_diel(
        axes_assignment_crys_to_diel::NTuple{3,Symbol};
        rotate_about::Union{Symbol,Nothing}=nothing,
        phi::Angle=0.0u"°"
    )

Returns the rotation matrix from crystallophysical (x,y,z) to dielectric (X,Y,Z) coordinates,
by optionally rotating around a dielectric axis after reassigning the axes.
- `axes_assignment_crys_to_diel`: 3-Tuple of Symbols (:X, :Y, :Z) indicating the assignment of crystallophysical to dielectric coordinates
- `rotate_about`: Symbol (:X, :Y, or :Z) indicating which dielectric axis to rotate around after assignment
- `phi`: rotation angle (applied after assignment)
"""
function rot_mat_crys_to_diel(
    axes_assignment_crys_to_diel::NTuple{3,Symbol};
    rotate_about::Union{Symbol,Nothing}=nothing,
    phi::Angle=0.0u"°"
)
    # Step 1: Assign crystallophysical axes (:x, :y, :z) to dielectric axes (:X, :Y, :Z)
    basis_crys = [
        SVector(1.0, 0.0, 0.0),  # x
        SVector(0.0, 1.0, 0.0),  # y
        SVector(0.0, 0.0, 1.0),  # z
    ]

    dielectric_axes = (:X, :Y, :Z)
    rot_mat = zeros(3, 3)
    for i in eachindex(dielectric_axes)
        crys_axis = axes_assignment_crys_to_diel[i]
        if crys_axis == :X
            idx = 1
        elseif crys_axis == :Y
            idx = 2
        elseif crys_axis == :Z
            idx = 3
        else
            error("Invalid crystallographic axis symbol: $crys_axis")
        end
        rot_mat[:, i] .= basis_crys[idx]
    end

    # Step 2: Apply rotation around dielectric axis after mapping
    if rotate_about !== nothing && phi != 0.0
        s, c = sincos(phi)
        if rotate_about == :X
            Rϕ = [
                1.0 0.0 0.0;
                0.0 c s;
                0.0 -s c
            ]
        elseif rotate_about == :Y
            Rϕ = [
                c 0.0 -s;
                0.0 1.0 0.0;
                s 0.0 c
            ]
        elseif rotate_about == :Z
            Rϕ = [
                c s 0.0;
                -s c 0.0;
                0.0 0.0 1.0
            ]
        else
            error("Invalid rotate_about axis: $rotate_about")
        end

        rot_mat = Rϕ * rot_mat
    end

    return rot_mat
end

function plot_axes_assignment_crys_to_diel(
    axes_assignment_crys_to_diel::NTuple{3,Symbol};
    rotate_about::Union{Symbol,Nothing}=nothing,
    phi::Angle=0.0u"°"
)
    R = rot_mat_crys_to_diel(axes_assignment_crys_to_diel; rotate_about, phi)

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, title="Crystallophysical (x, y, z) to dielectric (X, Y, Z) coordinates")
    hidedecorations!(ax)
    fac = 0.7

    # Original crystallophysical axes
    x = [1.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [0.0, 0.0, 1.0]
    crys_basis = [x, y, z]
    crys_labels = ("x", "y", "z")

    # Rotated dielectric axes
    diel_labels = ["X", "Y", "Z"]
    rot_axes = [R \ x, R \ y, R \ z]

    # Plot original axes (dashed blue)
    for (v, lbl) in zip(crys_basis, crys_labels)
        lines!(ax, [Point3f(0), Point3f(v)], color=:blue, linewidth=2, linestyle=:dash)
        text!(ax, v...; text=lbl, color=:blue)
    end

    # Plot rotated axes (solid red)
    for (v, lbl) in zip(rot_axes, diel_labels)
        lines!(ax, [Point3f(0), Point3f(fac * v)], color=:red, linewidth=2)
        text!(ax, (fac * v)...; text=lbl, color=:red)
    end

    return fig
end


function calc_d_XYZ_full(
    point_group::String,
    rot_mat::AbstractMatrix{<:Number}=I(3),
    use_kleinman::Bool=true;
    components_abc...
)
    comps = [c[1] for c in components_abc]
    vals = [c[2] for c in components_abc]
    sg = use_kleinman ? SYMMETRY_GROUPS_KLEINMAN : SYMMETRY_GROUPS

    symmetry = get(sg, point_group, nothing)
    symmetry !== nothing || error("Point group '$point_group' not defined.")

    # Iterate through each symmetry group to ensure consistency
    # Symmetry goups and given components are assumed to be specified in abc coordinates
    d_abc = zeros(3, 6) * u"pm/V"
    for group in symmetry
        group_comps, group_signs = group
        given_group_comps = intersect(Set(group_comps), Set(comps))

        if length(given_group_comps) > 1
            @error "Exactly one of the following components must be specified: $(group_comps)"
        elseif length(given_group_comps) < 1
            if length(group_comps) == 1
                @error "Component '$(group_comps[1])' must be specified."
            else
                @error "One of these components must be specified: $(group_comps)"
            end
        end

        # Correct signs within the current symmetry group
        given_group_comp = collect(given_group_comps)[1]
        idx = findall(c -> (c === given_group_comp), group_comps)[1]
        group_signs_switched = group_signs[idx] == 1 ? group_signs : -group_signs

        # Set all group components based on the given component with the correct sign relations
        c_idx = findall(c -> (c === given_group_comp), comps)[1]
        for (gc, gs) in zip(group_comps, group_signs_switched)
            d_abc[tensor_indices(gc)...] = vals[c_idx] * gs
        end
    end

    # Rotate d tensor from the coordinates system of the given d tensor components to XYZ coordinates based on the provided assigment
    d_abc_full = expand_d_contract(d_abc)
    d_XYZ_full = rotate_tensor3(d_abc_full, rot_mat)
    return d_XYZ_full
end

function calc_d_eff(
    cr::NonlinearCrystal,
    E_dir_r1::AbstractVector{<:Number},
    E_dir_r2::AbstractVector{<:Number},
    E_dir_b::AbstractVector{<:Number};
    lambda_rrb::Union{NTuple{3,Length},Nothing}=nothing,
    temp::Temperature=default_temp(cr),
)
    # Use Miller scaling if all lambdas are given and Miller's delta has been calculated during crystal initialization
    if isnothing(cr.miller_delta) || isnothing(lambda_rrb)
        d_XYZ_full = cr.d_XYZ_ref_full
    else
        d_XYZ_full = miller_rescale(cr, lambda_rrb; temp)
    end

    @tullio d_eff := E_dir_b[i] * d_XYZ_full[i, j, k] * E_dir_r1[j] * E_dir_r2[k]
    return d_eff
end

function miller_rescale(
    cr::NonlinearCrystal,
    lambda_rrb::Union{NTuple{3,Length},Nothing};
    temp::Temperature=default_temp(cr),
)
    d_XYZ_full = miller_rescale(
        cr.miller_delta,
        cr.n_X_principal,
        cr.n_Y_principal,
        cr.n_Z_principal,
        temp;
        lambda_r1=lambda_rrb[1],
        lambda_r2=lambda_rrb[2],
        lambda_b=lambda_rrb[3],
    )
    return d_XYZ_full
end

function miller_rescale(
    miller_delta::AbstractArray{<:Number,3},
    n_X_principal::RefractiveIndex,
    n_Y_principal::RefractiveIndex,
    n_Z_principal::RefractiveIndex,
    temp::Temperature;
    lambda_r1::Union{Length,Nothing}=nothing,
    lambda_r2::Union{Length,Nothing}=nothing,
    lambda_b::Union{Length,Nothing}=nothing,
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    # Compute refractive indices along XYZ at given wavelengths and temperature
    # Compute χ1 = n² - 1 values at reference measurement wavelengths
    n_funs = (n_X_principal, n_Y_principal, n_Z_principal)
    n_r1 = (f(lambda_r1, temp) for f in n_funs)
    n_r2 = (f(lambda_r2, temp) for f in n_funs)
    n_b = (f(lambda_b, temp) for f in n_funs)

    χ1_r1 = n_r1 .^ 2 .- 1
    χ1_r2 = n_r2 .^ 2 .- 1
    χ1_b = n_b .^ 2 .- 1

    # Miller’s Rule:
    # χ⁽²⁾ᵢⱼₖ = Δᵢⱼₖ · χ⁽¹⁾ᵢᵢ(ω_b) · χ⁽¹⁾ⱼⱼ(ω_r₁) · χ⁽¹⁾ₖₖ(ω_r₂)
    # where χ⁽¹⁾ₐₐ(ω) = nₐ²(ω) − 1
    @tullio d_XYZ_full[i, j, k] := miller_delta[i, j, k] * (χ1_b[i] * χ1_r1[j] * χ1_r2[k])
    return d_XYZ_full
end

function calc_miller_delta(
    d_ref_XYZ_full::AbstractArray{<:Number,3},
    n_X_principal::RefractiveIndex,
    n_Y_principal::RefractiveIndex,
    n_Z_principal::RefractiveIndex,
    temp_ref::Temperature;
    lambda_r1::Union{Length,Nothing}=nothing,
    lambda_r2::Union{Length,Nothing}=nothing,
    lambda_b::Union{Length,Nothing}=nothing,
)
    lambda_r1, lambda_r2, lambda_b = pm_wavelengths(; lambda_r1, lambda_r2, lambda_b)

    n_funs = (n_X_principal, n_Y_principal, n_Z_principal)
    n_r1 = (f(lambda_r1, temp_ref) for f in n_funs)
    n_r2 = (f(lambda_r2, temp_ref) for f in n_funs)
    n_b = (f(lambda_b, temp_ref) for f in n_funs)

    χ1_r1 = n_r1 .^ 2 .- 1
    χ1_r2 = n_r2 .^ 2 .- 1
    χ1_b = n_b .^ 2 .- 1

    # Miller’s Rule:
    # Δᵢⱼₖ = χ⁽²⁾ᵢⱼₖ / (χ⁽¹⁾ᵢᵢ(ω_b) · χ⁽¹⁾ⱼⱼ(ω_r₁) · χ⁽¹⁾ₖₖ(ω_r₂))
    # χ⁽¹⁾ₐₐ(ω) = nₐ²(ω) − 1
    @tullio miller_delta[i, j, k] := d_ref_XYZ_full[i, j, k] / (χ1_b[i] * χ1_r1[j] * χ1_r2[k])

    return SArray{Tuple{3,3,3},typeof(1.0u"m/V")}(miller_delta)
end

function calc_miller_delta(
    d_ref_XYZ_full::AbstractArray{<:Number,3},
    n_o_principal::RefractiveIndex,
    n_e_principal::RefractiveIndex,
    temp_ref::Temperature;
    lambda_r1::Union{Length,Nothing}=nothing,
    lambda_r2::Union{Length,Nothing}=nothing,
    lambda_b::Union{Length,Nothing}=nothing,
)
    n_X_principal = n_o_principal
    n_Y_principal = n_o_principal
    n_Z_principal = n_e_principal

    return calc_miller_delta(
        d_ref_XYZ_full,
        n_X_principal,
        n_Y_principal,
        n_Z_principal,
        temp_ref;
        lambda_r1,
        lambda_r2,
        lambda_b,
    )
end

function plot_miller_scaling_coeffs_shg(
    cr::NonlinearCrystal;
    temp::Temperature=default_temp(cr),
    size::NTuple{2,Int}=(800, 600),
)
    lambda_b_start = valid_lambda_range(cr)[1]
    lambda_b_end = valid_lambda_range(cr)[2] / 2
    all_lambda_b = LinRange(lambda_b_start:1.0u"nm":lambda_b_end)
    get_all_d = l -> abs.(miller_rescale(cr, (2 * l, 2 * l, l); temp))
    all_d = [contract_d_tensor(get_all_d(l)) for l in all_lambda_b]

    f = Figure(; size)
    uc1 = Makie.UnitfulConversion(u"µm"; units_in_label=true)
    uc2 = Makie.UnitfulConversion(u"pm/V"; units_in_label=true)
    ax = Axis(
        f[1, 1],
        xlabel="λ_b",
        ylabel="d_ij",
        dim1_conversion=uc1,
        dim2_conversion=uc2,
    )
    for i in CartesianIndices(all_d[1])
        all_d[1][i] == 0.0u"pm/V" && continue
        lab = (plot, index, position) -> "λ_b = $(round(u"µm", position[1] * u"µm"; digits=3))\nd$(i[1])$(i[2]) = $(round(u"pm/V", position[2] * u"pm/V"; digits=3))"
        lines!(ax, all_lambda_b, [d[i] for d in all_d], label="d$(i[1])$(i[2])", inspectable=true, inspector_label=lab)
    end
    DataInspector(ax)
    Legend(f[1, 2], ax)
    return f
end