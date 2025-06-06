module NonlinearCrystals

import Printf: @printf, @sprintf
using StaticArrays
using LinearAlgebra 
using Tullio 
using Unitful 
using Roots 
using GLMakie 
using GeometryBasics

import ForwardDiff
import PhysicalConstants.CODATA2022: c_0, ε_0

import Unitful: Temperature, Length, Frequency
@derived_dimension Angle  Unitful.𝐋^0  true # Used to force explicit angle units

# Color defaults
const COL_COORDS = to_colormap(:Greys)[7]
const COL_CONTOUR = to_colormap(:Set1_9)[1]
const COL_R1 = to_colormap(:vik10)[9]
const COL_R2 = to_colormap(:vik10)[8]
const COL_B = to_colormap(:vik10)[3]
const COLORMAP_HEATMAP = :vik

# Includes
include("utils.jl")
include("refractive_index.jl")
include("crystals.jl")
include("d_tensor.jl")
include("crystal_symmetry.jl")
include("phasematch.jl")
include("plot_phasematch.jl")

# Include all crystals
include("crystal_data/bbo.jl")
include("crystal_data/bibo.jl")
include("crystal_data/cga.jl")
include("crystal_data/gs.jl")
include("crystal_data/kdp.jl")
include("crystal_data/ktp_h.jl")
include("crystal_data/ktp_f.jl")
include("crystal_data/lbo.jl")
include("crystal_data/lnb_c.jl")
# include("crystal_data/lnb_m.jl") # TODO: Add when data is ready
include("crystal_data/lnb_s.jl")

# Global storage arrays
selected_pms = PhaseMatch[]
export selected_pms

end
