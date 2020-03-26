using RobotZoo
using Test
using StaticArrays
import RobotZoo: KnotPoint, dynamics, jacobian!, discrete_dynamics, discrete_jacobian!, RK3

@testset "allocations" begin
    include("models_test.jl")
end
