module RobotZoo

using StaticArrays
using Parameters
using RobotDynamics
using Rotations
using LinearAlgebra

using RobotDynamics: ContinuousDynamics, RigidBody, LieGroupModel, ContinuousDynamics
import RobotDynamics: dynamics!, dynamics, forces, moments, wrenches, inertia, inertia_inv, orientation
import RobotDynamics: state_dim, control_dim

include("acrobot.jl")
include("car.jl")
include("cartpole.jl")
include("double_integrator.jl")
include("pendulum.jl")
include("quadrotor.jl")
include("planar_quad.jl")
include("yak_plane.jl")
include("satellite.jl")
include("freebody.jl")
include("rocket.jl")

include("LinearModels.jl")

end # module
