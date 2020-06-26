"""
    Pendulum{R}

The most simple canonical nonlinear system, consisting of a simple pendulum controlled by
torque at the base. Has 2 state and 1 control.

# Constructor
    Pendulum(; kwargs...)

with keyword arguments
* `mass` - mass of the pendulum, in kg (default = 1.0)
* `length` - length of the pendulum, in m (default = 0.5)
* `b` - friction coefficient (default=0.1)
* `I` - inertia, in kg⋅m² (default=0.25)
* `g` - gravity, in kg/m² (default=9.81)
"""
@with_kw mutable struct Pendulum{T} <: AbstractModel
    mass::T = 1.
    length::T = 0.5
    b::T = 0.1
    lc::T = 0.5
    I::T = 0.25
    g::T = 9.81
end

RobotDynamics.state_dim(::Pendulum) = 2
RobotDynamics.control_dim(::Pendulum) = 1

function dynamics(p::Pendulum, x, u)
    m = p.mass * p.lc * p.lc
    @SVector [x[2],
              u[1]/m - p.g*sin(x[1])/p.lc - p.b*x[2]/m]
end
@inline Base.position(::Pendulum, x) = @SVector zeros(3)
orientation(::Pendulum, x) = expm((pi-x[1])*@SVector [1,0,0.])
