"""
    DubinsCar

A simple unicyle model with state vector  `[x,y,θ]` and control vector `[v,ω]`.

# Constructor
    DubinsCar(; [radius=0.175])

where `radius` is the radius of the car, used only for plotting and obstacle avoidance constraints.
"""
DubinsCar
@autodiff struct DubinsCar <: ContinuousDynamics
    radius::Float64  # [m] radius of a Roomba
    function DubinsCar(;radius=0.175)
        new(radius)
    end
end
RobotDynamics.state_dim(::DubinsCar) = 3
RobotDynamics.control_dim(::DubinsCar) = 2

function dynamics(::DubinsCar,x,u)
    ẋ = @SVector [u[1]*cos(x[3]),
                  u[1]*sin(x[3]),
                  u[2]]
end
function dynamics!(::DubinsCar,ẋ,x,u)
    ẋ[1] = u[1]*cos(x[3])
    ẋ[2] = u[1]*sin(x[3])
    ẋ[3] = u[2]
    return nothing
end

function jacobian!(::DubinsCar, J,xdot,x,u)
    sθ,cθ = sincos(x[3])
    J .= 0
    J[1,3] = -u[1]*sθ
    J[1,4] = cθ 
    J[2,3] = u[1]*cθ
    J[2,4] = sθ
    J[3,5] = 1.0
    return nothing
end

Base.position(::DubinsCar, x) = @SVector [x[1], x[2], 0.0]
orientation(::DubinsCar, x) = expm(x[3]*@SVector [0,0,1.])

"""
    BicycleModel <: ContinuousDynamics

Kinematic model of a car with front-wheel steering. Assumes no skid and models wheels on the same
axel as a single wheel. The state is defined as `[x,y,θ,δ]` where `x` and `y` are the position,
`θ` is the orientation (yaw angle), and `δ` is the steering angle. The controls are `[v,ϕ]` where 
`v` is the forward velocity and `ϕ` is the steering angle rate.

The reference location on the car can be specified using the `ref` field, which is either at 
the center of the rear wheel `:rear`, the front wheel `:front`, or at the center of gravity `:cg` (default).

# Constructor
    BicycleModel(;kwargs...)

with keyword arguments:
* `ref`: location of the reference position on the vehicle. Valid values: `[:rear, :front, :cg]`
* `L`: distance between the front and rear wheels
* `lr`: distance from the center of gravity to the center of the rear wheel
"""
@autodiff struct BicycleModel <: ContinuousDynamics
    ref::Symbol
    L::Float64    # distance between wheels  (m)
    lr::Float64   # distance to rear wheels  (m)
    function BicycleModel(;ref::Symbol=:cg, L::Real=2.7, lr::Real=1.5)
        @assert ref ∈ (:cg, :front, :rear)
        @assert L > 0 "Vehicle length ($L) must be greater than 0"
        @assert L > lr "Vehicle length ($L) must be greater than the distance to the rear wheel ($lr)"
        new(ref, L, lr)
    end
end

fun = :foo
args = [:(model::BicycleModel), :x, :u]

body = quote
    v = u[1]  # longitudinal velocity (cmd)
    ϕ = u[2]  # steering angle rate   (cmd)
    θ = x[3]  # yaw
    δ = x[4]  # steering angle
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
end
@eval function dynamics(model::BicycleModel, x, u)
    $body
    return SA[ẋ, ẏ, ω, ϕ]
end
@eval function dynamics!(model::BicycleModel, xdot, x, u)
    $body
    xdot[1] = ẋ
    xdot[2] = ẏ
    xdot[3] = ω
    xdot[4] = ϕ
    return nothing
end

state_dim(::BicycleModel) = 4
control_dim(::BicycleModel) = 2


# using Plots
# function Plots.plot(model::DubinsCar, Z::Traj; kwargs...)
#     x = [z.z[1] for z in Z]
#     y = [z.z[2] for z in Z]
#     Plots.plot(x, y, xlabel="x position", ylabel="y position", markershape=:circle;
#         kwargs...)
# end
# function Plots.plot!(model::DubinsCar, Z::Traj; kwargs...)
#     x = [z.z[1] for z in Z]
#     y = [z.z[2] for z in Z]
#     Plots.plot!(x, y, xlabel="x position", ylabel="y position", markershape=:circle;
#         kwargs...)
# end
