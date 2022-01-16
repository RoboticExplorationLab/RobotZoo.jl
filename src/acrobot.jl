
"""
    Acrobot{T}

A double-pendulum with actuation only at the elbow joint. It has 4 states and 1 control.

# Constructor
    Acrobot(; kwargs...)

with keyword arguments:
* `l` - `SVector{2,T}` of link lengths
* `m` - `SVector{2,T}` of link masses
* `J` - `SVector{2,T}` of link inertias
"""
@autodiff struct Acrobot{T} <: ContinuousDynamics
    l::SVector{2,T}
    m::SVector{2,T}
    J::SVector{2,T}
end

function Acrobot(;l=SA[1.,1.],m=SA[1.,1.],J=SA[1.,1.])
    T = promote_type(eltype(l), eltype(m), eltype(J)) 
    Acrobot{T}(SA[l[1],l[2]], SA[m[1],m[2]], SA[J[1],J[2]])
end

function Acrobot(l,m,J)
    T = eltype(l)
    Acrobot(SA[l[1],l[2]], SA[m[1],m[2]], SA[J[1],J[2]])
end

function dynamics(model::Acrobot, x, u)
    g = 9.81
    m1,m2 = model.m
    l1,l2 = model.l
    J1,J2 = model.J
    θ1,    θ2    = x[1], x[2]
    θ1dot, θ2dot = x[3], x[4]
    s1,c1 = sincos(θ1)
    s2,c2 = sincos(θ2)
    c12 = cos(θ1 + θ2)

    # mass matrix
    m11 = m1*l1^2 + J1 + m2*(l1^2 + l2^2 + 2*l1*l2*c2) + J2
    m12 = m2*(l2^2 + l1*l2*c2 + J2)
    m22 = l2^2*m2 + J2
    M = @SMatrix [m11 m12; m12 m22]

    # bias term
    tmp = l1*l2*m2*s2
    b1 = -(2 * θ1dot * θ2dot + θ2dot^2)*tmp
    b2 = tmp * θ1dot^2
    B = @SVector [b1, b2]

    # friction
    c = 1.0
    C = @SVector [c*θ1dot, c*θ2dot]

    # gravity term
    g1 = ((m1 + m2)*l2*c1 + m2*l2*c12) * g
    g2 = m2*l2*c12*g
    G = @SVector [g1, g2]

    # equations of motion
    τ = @SVector [0, u[1]]
    θddot = M\(τ - B - G - C)
    return @SVector [θ1dot, θ2dot, θddot[1], θddot[2]]
end

function dynamics!(model, xdot, x, u)
    xdot .= dynamics(model, x, u)
end

RobotDynamics.state_dim(::Acrobot) = 4
RobotDynamics.control_dim(::Acrobot) = 1
