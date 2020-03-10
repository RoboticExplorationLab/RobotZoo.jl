module RobotZoo

using StaticArrays
using Parameters
using RobotDynamics
using DifferentialRotations
using LinearAlgebra

import Dynamics: dynamics, forces, moments, wrenches, mass_matrix, inertia, inertia_inv, orientation
import Dynamics: state_dim, control_dim

include("acrobot.jl")
include("car.jl")
include("cartpole.jl")
include("double_integrator.jl")
include("pendulum.jl")
include("quadrotor.jl")
include("yak_plane.jl")

end # module
