export crystal_system, plot_symmetry

# Mapping from point group to crystal system
const POINT_GROUP_MAP = Dict(
    # Triclinic
    "1"     => "Triclinic",
    "-1"    => "Triclinic",

    # Monoclinic
    "2"     => "Monoclinic",
    "2/m"   => "Monoclinic",
    "m"     => "Monoclinic",

    # Orthorhombic
    "222"   => "Orthorhombic",
    "mm2"   => "Orthorhombic",
    "mmm"   => "Orthorhombic",

    # Tetragonal
    "4"     => "Tetragonal",
    "-4"    => "Tetragonal",
    "4/m"   => "Tetragonal",
    "422"   => "Tetragonal",
    "4mm"   => "Tetragonal",
    "-42m"  => "Tetragonal",
    "4/mmm" => "Tetragonal",

    # Trigonal (Rhombohedral)
    "3"     => "Trigonal",
    "-3"    => "Trigonal",
    "32"   => "Trigonal",
    "3m"    => "Trigonal",
    "-3m"   => "Trigonal",

    # Hexagonal
    "6"     => "Hexagonal",
    "-6"    => "Hexagonal",
    "6/m"   => "Hexagonal",
    "622"   => "Hexagonal",
    "6mm"   => "Hexagonal",
    "-6m2"  => "Hexagonal",
    "6/mmm" => "Hexagonal",

    # Cubic
    "23"    => "Cubic",
    "m-3"    => "Cubic",
    "432"   => "Cubic",
    "-43m"  => "Cubic",
    "m-3m"   => "Cubic"
)

"""
    crystal_system(point_group::AbstractString)

Returns the name of the crystal system corresponding to the given crystallographic point group (e.g., `"mm2" → "Orthorhombic"`).
Throws an error if the point group is unknown. Based on conventional assignments from nonlinear optics literature.
"""
function crystal_system(point_group::AbstractString) 
    if haskey(POINT_GROUP_MAP, point_group)
        return POINT_GROUP_MAP[point_group]
    else
        @error "Point group '$(point_group)' is unknown."
    end
end

crystal_system(cr::NonlinearCrystal) = crystal_system(cr.metadata[:point_group])

# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#zernike06applied
const SYMMETRY_GROUPS_KLEINMAN = Dict(
    # Pairs of related coefficients based on Kleinman's symmetry conditions and sign relations within each component group 
    "1" => [
        ([:d11], [1]),
        ([:d21, :d16], [1, 1]),
        ([:d31, :d15], [1, 1]),
        ([:d12, :d26], [1, 1]),
        ([:d22], [1]),
        ([:d32, :d24], [1, 1]),
        ([:d13, :d35], [1, 1]),
        ([:d23, :d34], [1, 1]),
        ([:d33], [1]),
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "2" => [
        ([:d21, :d16], [1, 1]),
        ([:d22], [1]),
        ([:d13], [1]),
        ([:d23, :d34], [1, 1]),
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "m" => [
        ([:d11], [1]),
        ([:d31, :d15], [1, 1]),
        ([:d12, :d26], [1, 1]),
        ([:d32, :d24], [1, 1]),
        ([:d13, :d35], [1, 1]),
    ],
    "222" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "mm2" => [
        ([:d31, :d15], [1, 1]),
        ([:d32, :d24], [1, 1]),
        ([:d33], [1]),
    ],
    "3" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
        ([:d21, :d22, :d16], [-1, 1, -1]),
        ([:d31, :d32, :d24, :d15], [1, 1, 1, 1]),
        ([:d33], [1]),
    ],
    "3m" => [
        ([:d21, :d22, :d16], [-1, 1, -1]),
        ([:d31, :d32, :d24, :d15], [1, 1, 1, 1]),
        (([:d33], [1]))
    ],
    "-6" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
        ([:d21, :d22, :d16], [-1, 1, -1]),
    ],
    "-6m2" => [
        ([:d21, :d22, :d16], [-1, 1, -1]),
    ],
    "6" => [
        ([:d31, :d32, :d24, :d15], [1, -1, 1, 1]),
        ([:d33], [1]),
    ],
    "4" => [
        ([:d31, :d32, :d24, :d15], [1, -1, 1, 1]),
        ([:d33], [1]),
    ],
    "6mm" => [
        ([:d31, :d32, :d24, :d15], [1, -1, 1, 1]),
        ([:d33], [1]),
    ],
    "4mm" => [
        ([:d31, :d32, :d24, :d15], [1, -1, 1, 1]),
        ([:d33], [1]),
    ],
    "622" => [([], [])],
    "422" => [([], [])],
    "-4" => [
        ([:d31, :d32, :d24, :d15], [1, -1, -1, 1]),
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "32" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
    ],
    "-42m" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "-43m" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "23" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "432" => [([], [])],
)

# Pairs of related coefficients due to the point group and sign relations within each component group 
# Reference: https://martinkosch.github.io/NonlinearCrystals.jl/dev/bibliography/#zernike06applied
const SYMMETRY_GROUPS = Dict(
    "1" => [
        ([:d11], [1]),
        ([:d21], [1]),
        ([:d31], [1]),
        ([:d12], [1]),
        ([:d22], [1]),
        ([:d32], [1]),
        ([:d13], [1]),
        ([:d23], [1]),
        ([:d33], [1]),
        ([:d14], [1]),
        ([:d24], [1]),
        ([:d34], [1]),
        ([:d15], [1]),
        ([:d25], [1]),
        ([:d35], [1]),
        ([:d16], [1]),
        ([:d26], [1]),
        ([:d36], [1]),
    ],
    "2" => [
        ([:d21], [1]),
        ([:d16], [1]),
        ([:d22], [1]),
        ([:d13], [1]),
        ([:d23], [1]),
        ([:d34], [1]),
        ([:d14], [1]),
        ([:d25], [1]),
        ([:d36], [1]),
    ],
    "m" => [
        ([:d11], [1]),
        ([:d31], [1]),
        ([:d15], [1]),
        ([:d12], [1]),
        ([:d26], [1]),
        ([:d32], [1]),
        ([:d24], [1]),
        ([:d13], [1]),
        ([:d35], [1]),
    ],
    "222" => [
        ([:d14], [1]),
        ([:d25], [1]),
        ([:d36], [1])
    ],
    "mm2" => [
        ([:d31], [1]),
        ([:d32], [1]),
        ([:d33], [1]),
        ([:d24], [1]),
        ([:d15], [1]),
    ],
    "3" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
        ([:d21, :d22, :d16], [-1, 1, -1]),
        ([:d31, :d32], [1, 1]),
        ([:d24, :d15], [1, 1]),
        ([:d33], [1]),
        ([:d14, :d25], [1, -1]),
    ],
    "3m" => [
        ([:d21, :d22, :d16], [-1, 1, -1]),
        ([:d31, :d32], [1, 1]),
        ([:d24, :d15], [1, 1]),
        (([:d33], [1]))
    ],
    "-6" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
        ([:d21, :d22, :d16], [-1, 1, -1]),
    ],
    "-6m2" => [
        ([:d21, :d22, :d16], [-1, 1, -1]),
    ],
    "6" => [
        ([:d31, :d32], [1, -1]),
        ([:d24, :d15], [1, 1]),
        ([:d33], [1]),
        ([:d14, :d25], [1, -1]),
    ],
    "4" => [
        ([:d31, :d32], [1, -1]),
        ([:d24, :d15], [1, 1]),
        ([:d33], [1]),
        ([:d14, :d25], [1, -1]),
    ],
    "6mm" => [
        ([:d31, :d32], [1, -1]),
        ([:d24, :d15], [1, 1]),
        ([:d33], [1]),
    ],
    "4mm" => [
        ([:d31, :d32], [1, -1]),
        ([:d24, :d15], [1, 1]),
        ([:d33], [1]),
    ],
    "622" => [
        ([:d14, :d25], [1, -1]),
    ],
    "422" => [
        ([:d14, :d25], [1, -1]),
    ],
    "-4" => [
        ([:d31, :d32], [1, -1]),
        ([:d24, :d15], [-1, 1]),
        ([:d14, :d25], [1, 1]),
        ([:d36], [1]),
    ],
    "32" => [
        ([:d11, :d12, :d26], [1, -1, -1]),
        ([:d14, :d25], [1, -1]),
    ],
    "-42m" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "-43m" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "23" => [
        ([:d14, :d25, :d36], [1, 1, 1]),
    ],
    "432" => [([], [])],
)

function find_zero_under_kleinman(pg::String)
    # Gather all coefficients from normal symmetry
    normal_coefs = Set{Symbol}()
    for (coefs, _) in get(SYMMETRY_GROUPS, pg, [])
        union!(normal_coefs, coefs)
    end

    # Gather all coefficients from Kleinman symmetry
    kleinman_coefs = Set{Symbol}()
    for (coefs, _) in get(SYMMETRY_GROUPS_KLEINMAN, pg, [])
        union!(kleinman_coefs, coefs)
    end

    # Anything in normal but missing from Kleinman: zero under Kleinman
    return setdiff(normal_coefs, kleinman_coefs)
end

"""
    plot_symmetry(point_group::AbstractString)

Visualizes the symmetry relations between nonlinear tensor components `d_ij` for a given point group.

Black lines indicate symmetry relations that follow from the basic point group; dashed lines show those that appear only under Kleinman symmetry assumptions.
Black circles represent same-sign components; white circles represent opposite signs. Square markers indicate coefficients that become zero under Kleinman conditions.
"""
function plot_symmetry(point_group::AbstractString)
    # Coordinates for each d_ij in a 6×3 grid
    coords = Dict(
        :d11 => (1, 3), :d12 => (2, 3), :d13 => (3, 3), :d14 => (4, 3), :d15 => (5, 3), :d16 => (6, 3),
        :d21 => (1, 2), :d22 => (2, 2), :d23 => (3, 2), :d24 => (4, 2), :d25 => (5, 2), :d26 => (6, 2),
        :d31 => (1, 1), :d32 => (2, 1), :d33 => (3, 1), :d34 => (4, 1), :d35 => (5, 1), :d36 => (6, 1)
    )

    fig = Figure(size=(600, 300))
    ax = Axis(fig[1, 1], title="Class: $point_group")

    # Mark all points as small black dots initially
    for (comp, (x, y)) in coords
        scatter!(ax, [x], [y], color=:black, markersize=6)
    end

    # Compute the set of coefficients that vanish under Kleinman
    zero_k_set = find_zero_under_kleinman(point_group)

    # Helper: draw connections
    function draw_connections(groups, linestyle=:solid; is_kleinman=false)
        for (group, signs) in groups
            xy = [coords[c] for c in group if haskey(coords, c)]
            xvals, yvals = getindex.(xy, 1), getindex.(xy, 2)

            # Draw lines
            for i in 1:length(xvals)-1
                lines!(ax, [xvals[i], xvals[i+1]], [yvals[i], yvals[i+1]],
                    color=:black, linestyle=linestyle, linewidth=2)
            end

            # Draw each point
            for (comp, sgn) in zip(group, signs)
                # If coefficient is zero under Kleinman, show a square
                shape = comp in zero_k_set ? :rect : :circle
                fillcolor = (sgn == 1) ? :black : :white
                size = 16
                (cx, cy) = coords[comp]
                scatter!(ax, [cx], [cy],
                    markersize=size, color=fillcolor, marker=shape,
                    strokewidth=1.5, strokecolor=:black)
            end
        end
    end

    # 1) Always‐valid symmetries (solid lines)
    main_groups = get(SYMMETRY_GROUPS, point_group, [])
    draw_connections(main_groups, :solid; is_kleinman=false)

    # 2) Kleinman‐only symmetries (dashed lines)
    all_kgroups = get(SYMMETRY_GROUPS_KLEINMAN, point_group, [])
    # remove any that appear identically in main_groups
    # so only truly "extra" Kleinman sets get drawn as dashed
    unique_kleinman_groups = setdiff(all_kgroups, main_groups)
    draw_connections(unique_kleinman_groups, :dash; is_kleinman=true)

    sym_ls = LineElement(color=:black, linestyle=:solid, linewidth=2)
    sym_kleinman = LineElement(color=:black, linestyle=:dash, linewidth=2)
    dot_plus = MarkerElement(color=:black, marker=:circle, markersize=10)
    dot_minus = MarkerElement(color=:white, marker=:circle, strokewidth=1.5, strokecolor=:black, markersize=10)
    square_zero_k_plus = MarkerElement(color=:black, marker=:rect, strokecolor=:black, markersize=10)
    square_zero_k_minus = MarkerElement(color=:white, marker=:rect, strokewidth=1.5, strokecolor=:black, markersize=10)

    Legend(fig[1, 2],
        [sym_ls, sym_kleinman, dot_plus, dot_minus, square_zero_k_plus, square_zero_k_minus],
        ["Always valid", "Kleinman only", "Same sign", "Opposite sign", "Zero under Kleinman, same sign", "Zero under Kleinman, opposite sign"],
        framevisible=false
    )

    hidexdecorations!(ax)
    hideydecorations!(ax)
    hidespines!(ax)

    fig
end

plot_symmetry(cr::NonlinearCrystal) = plot_symmetry(cr.metadata[:point_group])