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
    function Rover(;ref::Symbol=:rear, L::Real=2.7, lr::Real=1.5, accel_lim::Real=2.5, steering_limit::Real=0.1)
        @assert ref ∈ (:cg, :front, :rear)
        @assert L > 0 "Vehicle length ($L) must be greater than 0"
        @assert L > lr "Vehicle length ($L) must be greater than the distance to the rear wheel ($lr)"
        new(ref, L, lr, accel_lim, steering_limit)
    end

end

fun = :foo
args = [:(model::Rover), :x, :u]

body = quote
    a_in = u[1]  # longitudinal acceleration (cmd)
    δ_in = u[2]  # steering angle (cmd)     check order

    if (-model.accel_lim > a_in)
        a_in = -model.accel_lim
    elseif a_in > model.accel_lim
        a_in = model.accel_lim
    end

    if (-model.steering_limit > δ_in)
        δ_in = -model.steering_limit
    elseif δ_in > model.steering_limit
        δ_in = model.steering_limit
    end

    s = x[1]
    y = x[2]
    θ = x[3]  # yaw
    δ = x[4]  # steering angle
    v = x[5]
    a = x[6]

    # Dynamics
    λ₁ = 100
    λ₂ = 100
    c = cos(θ)
    s = sin(θ)

    ẋ = v*c
    ẏ = v*s
    ω = v*tan(δ_in)       # yaw rate
    ϕ = -λ₁*(δ - δ_in) # steering angle rate
    # ϕ = 0
    v̇ = a_in
    ȧ = -λ₂*(a - a_in)
    # ȧ = 0
end
@eval function dynamics(model::Rover, x, u)
    $body
    return SA[ẋ, ẏ, ω, ϕ, v̇, ȧ]
end
@eval function dynamics!(model::Rover, xdot, x, u)
    $body
    xdot[1] = ẋ
    xdot[2] = ẏ
    xdot[3] = ω
    xdot[4] = ϕ
    xdot[5] = v̇
    xdot[6] = ȧ
    return nothing
end

# using Rotations
# using GeometryBasics: HyperRectangle, Vec, Point, Mesh
# using MeshCat
# using CoordinateTransformations
# function visualize(filepath, xs)
#     vis = Visualizer()
#     render(vis)
#     anim = Animation()
#     rover_base_mesh = MeshFileGeometry(filepath + "/base.obj")
#     rover_wheel_left = MeshFileGeometry("./Rover Model/frontwheel.obj")
#     rover_wheel_right = MeshFileGeometry("./Rover Model/frontwheel.obj")
#     setobject!(vis[:rc][:ori][:rover], rover_base_mesh)
#     setobject!(vis[:rc][:ori][:rover][:wheelright][:steer], rover_wheel_right)
#     setobject!(vis[:rc][:ori][:rover][:wheelleft][:steer], rover_wheel_left)
#     settransform!(vis[:rc], Translation(0, -0.075, 0) ∘ LinearMap(RotX(π/2)))
#     settransform!(vis[:rc][:ori], LinearMap(RotY(π)))
#     settransform!(vis[:rc][:ori][:rover][:wheelleft], Translation(0.035, 0.06, 0.025) ∘ LinearMap(RotX(π)))
#     settransform!(vis[:rc][:ori][:rover][:wheelright], Translation(0.035, 0, .125))
# end

state_dim(::Rover) = 6
control_dim(::Rover) = 2