module RobotZoo

using StaticArrays
using Parameters
using RobotDynamics
using Rotations
using LinearAlgebra

import RobotDynamics: dynamics, forces, moments, wrenches, inertia, inertia_inv, orientation
import RobotDynamics: state_dim, control_dim

include("acrobot.jl")
include("car.jl")
include("cartpole.jl")
include("double_integrator.jl")
include("pendulum.jl")
include("quadrotor.jl")
include("yak_plane.jl")
include("satellite.jl")

end # module
