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
