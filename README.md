![CI](https://github.com/bjack205/RobotZoo.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/RoboticExplorationLab/RobotZoo.jl/branch/master/graph/badge.svg?token=KWcOu2TFpd)](https://codecov.io/gh/RoboticExplorationLab/RobotZoo.jl)

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
import RobotDynamics as RD

model = RobotZoo.Cartpole()
n,m = RD.dims(model)

# Generate random state and control vector
x,u = rand(model)
t = 0.0   # time (s)
dt = 0.1  # time step (s)
z = RD.KnotPoint(x,u,t,dt)

# Evaluate the continuous dynamics and Jacobian
ẋ = dynamics(model, x, u)
∇f = zeros(n, n + m)
RD.jacobian!(RD.StaticReturn(), RD.ForwardAD(), model, ∇f, ẋ, model, z)

# Evaluate the discrete dynamics and Jacobian
dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
x′ = RD.discrete_dynamics(dmodel, model, x, u, t, dt)
RD.discrete_jacobian!(RD.StaticReturn(), RD.ForwardAD(), dmodel, ∇f, x′, z)
```
