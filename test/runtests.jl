using NonlinearCrystals
using Test
using Aqua

@testset "NonlinearCrystals.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(NonlinearCrystals)
    end
    # Write your tests here.
end
