"""
Rover <: ContinuousDynamics

Kinematic model of a car with front-wheel steering. Assumes no skid and models wheels on the same
axel as a single wheel. The state is defined as `[s,y,θ,δ,v,a]` where `s` and `y` are the position,
`θ` is the orientation (yaw angle), `δ` is the steering angle, `v` is the longitudinal velocity, and
`a` is the longitudinal acceleration. The controls are `[a_in,δ_in]` where `a_in` is the forward 
acceleration and `δ_in` is the steering angle.

The reference location on the car can be specified using the `ref` field, which is either at 
the center of the rear wheel `:rear`, the front wheel `:front`, or at the center of gravity `:cg` (default).

# Constructor
Rover(;kwargs...)

with keyword arguments:
* `ref`: location of the reference position on the vehicle. Valid values: `[:rear, :front, :cg]`
* `L`: distance between the front and rear wheels
* `lr`: distance from the center of gravity to the center of the rear wheel
"""
Rover
@autodiff struct Rover <: ContinuousDynamics
    ref::Symbol
    L::Float64    # distance between wheels  (m)
    lr::Float64   # distance to rear wheels  (m)
    accel_lim::Float64
    steering_limit::Float64
    function Rover(;ref::Symbol=:rear, L::Real=2.7, lr::Real=1.5, accel_lim::Real=4.0, steering_limit::Real=.15)
        @assert ref ∈ (:cg, :front, :rear)
        @assert L > 0 "Vehicle length ($L) must be greater than 0"
        @assert L > lr "Vehicle length ($L) must be greater than the distance to the rear wheel ($lr)"
        new(ref, L, lr, accel_lim, steering_limit)
    end
end

fun = :foo
args = [:(model::Rover), :x, :u]

body = quote
    α = u[1]  # longitudinal acceleration (cmd)
    ϕ = u[2]  # steering angle rate (cmd)

    if (-model.accel_lim > α)
        α = -model.accel_lim
    elseif α > model.accel_lim
        α = model.accel_lim
    end

    if (-model.steering_limit > ϕ)
        ϕ = -model.steering_limit
    elseif ϕ > model.steering_limit
        ϕ = model.steering_limit
    end

    θ = x[3]  # yaw
    δ = x[4]  # steering angle
    v = x[5]

    # Dynamics
    # c = cos(θ)
    # s = sin(θ)
    
    if model.ref == :cg
        β = atan(model.lr * δ, model.L)
        s,c = sincos(θ + β)
        ω = v*cos(β)*tan(δ) / model.L
    elseif model.ref == :rear
        s,c = sincos(θ)
        ω = v*tan(δ) / model.L
    elseif model.ref == :front
        s,c = sincos(θ + δ)
        ω = v*sin(δ) / model.L
    end

    ẋ = v*c
    ẏ = v*s
    v̇ = α
end
@eval function dynamics(model::Rover, x, u)
    $body
    return SA[ẋ, ẏ, ω, ϕ, v̇]
end
@eval function dynamics!(model::Rover, xdot, x, u)
    $body
    xdot[1] = ẋ
    xdot[2] = ẏ
    xdot[3] = ω
    xdot[4] = ϕ
    xdot[5] = v̇
    return nothing
end

state_dim(::Rover) = 5
control_dim(::Rover) = 2