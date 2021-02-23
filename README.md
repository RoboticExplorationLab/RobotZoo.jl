![CI](https://github.com/bjack205/RobotZoo.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/bjack205/RobotZoo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/bjack205/RobotZoo.jl)

# RobotZoo.jl

Provides a handful of dynamics models of common robotic platforms, implemented using the 
[`RobotDynamics.jl`](https://github.com/RoboticExplorationLab/RobotDynamics.jl) package. 
The following models are currently implemented:
* Acrobot (`Acrobot`)
* Dubins Car (unicycle model) (`DubinsCar`)
* Kinematic Bicycle car model (`BicycleModel`)
* Cartpole (`Cartpole`)
* Double Integrator (`DoubleIntegrator`)
* Pendulum (`Pendulum`)
* Quadrotor (`Quadrotor`)
* Airplane (`YakPlane`)
* Satellite (`Satellite`)

Most models can be constructed using their default constructor, e.g. `model = RobotZoo.Cartpole()`.
To get more information on each model, refer to the documentation for each type, accessible via the command line:
```julia
?RobotZoo.Quadrotor
```

## Example Usage
```julia
using RobotZoo
using RobotDynamics

model RobotZoo.Cartpole()
n,m = size(model)

# Generate random state and control vector
x,u = rand(model)
dt = 0.1  # time step (s)
z = KnotPoint(x,u,dt)

# Evaluate the continuous dynamics and Jacobian
ẋ = dynamics(model, x, u)
∇f = RobotDynamics.DynamicsJacobian(model)
jacobian!(∇f, model, z)

# Evaluate the discrete dynamics and Jacobian
x′ = discrete_dynamics(RK3, model, z)
discrete_jacobian!(RK3, ∇f, model, z)
```
