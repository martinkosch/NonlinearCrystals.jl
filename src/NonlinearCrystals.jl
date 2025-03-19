module NonlinearCrystals

using StaticArrays
using LinearAlgebra 
using Unitful 
using Roots 
using GLMakie 
using GeometryBasics

import ForwardDiff
import PhysicalConstants.CODATA2022: c_0

include("utils.jl")
include("refractive_index.jl")
include("crystal_symmetry.jl")
include("crystals.jl")
include("phasematch.jl")

# Include all crystals
include("crystal_data/bbo.jl")
include("crystal_data/cga.jl")
include("crystal_data/ktp_h.jl")
include("crystal_data/ktp_f.jl")

end
