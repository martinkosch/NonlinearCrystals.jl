using NonlinearCrystals
using SafeTestsets
using Test
using Aqua

@safetestset "Aqua" begin
    include("aqua.jl")
end

@safetestset "Utils" begin
    include("utils.jl")
end


@safetestset "Refractive index" begin
    include("refractive_index.jl")
end


## Crystal data tests

@safetestset "Crystal test: BBO" begin
    include("crystal_data_tests/bbo.jl")
end

@safetestset "Crystal test: CGA" begin
    include("crystal_data_tests/cga.jl")
end

@safetestset "Crystal test: KTP_F" begin
    include("crystal_data_tests/ktp_f.jl")
end

@safetestset "Crystal test: KTP_H" begin
    include("crystal_data_tests/ktp_h.jl")
end

@safetestset "Crystal test: LBO" begin
    include("crystal_data_tests/lbo.jl")
end