using RobotZoo
using RobotDynamics
using Test
using StaticArrays
using BenchmarkTools
using Rotations

@testset "allocations" begin
    include("models_test.jl")
end

@testset "LinearModels" begin
    include("linear_models_test.jl")
end
