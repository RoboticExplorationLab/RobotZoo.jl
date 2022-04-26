using RobotZoo
using RobotDynamics
using Test
using StaticArrays
using BenchmarkTools
using Rotations
using Random

@testset "Allocations" begin
    include("models_test.jl")
end

@testset "LinearModels" begin
    include("linear_models_test.jl")
    include("linear_models_examples_test.jl")
end
