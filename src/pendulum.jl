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

@autodiff mutable struct Pendulum{T} <: ContinuousDynamics
    mass::T
    len::T
    b::T
    lc::T
    I::T
    g::T
end
function Pendulum(;mass=1., len=0.5, b=0.1, lc=0.5, I=0.25, g=9.81)
    T = eltype(promote(mass, len, b, lc, I, g))
    Pendulum{T}(mass, len, b, lc, I, g)
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
