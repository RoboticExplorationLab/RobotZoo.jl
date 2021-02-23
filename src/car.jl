"""
    DubinsCar

A simple unicyle model with state vector  `[x,y,θ]` and control vector `[v,ω]`.

# Constructor
    DubinsCar(; [radius=0.175])

where `radius` is the radius of the car, used only for plotting and obstacle avoidance constraints.
"""
@with_kw struct DubinsCar <: AbstractModel
    radius::Float64 = 0.175  # [m] radius of a Roomba
end
RobotDynamics.state_dim(::DubinsCar) = 3
RobotDynamics.control_dim(::DubinsCar) = 2

function dynamics(::DubinsCar,x,u)
    ẋ = @SVector [u[1]*cos(x[3]),
                  u[1]*sin(x[3]),
                  u[2]]
end

Base.position(::DubinsCar, x) = @SVector [x[1], x[2], 0.0]
orientation(::DubinsCar, x) = expm(x[3]*@SVector [0,0,1.])

"""
    BicycleModel <: AbstractModel

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
@with_kw struct BicycleModel <: AbstractModel
    ref::Symbol = :cg
    L::Float64  = 2.7   # distance between wheels  (m)
    lr::Float64 = 1.5   # distance to rear wheels  (m)
    function BicycleModel(ref::Symbol, L::Real, lr::Real)
        @assert ref ∈ (:cg, :front, :rear)
        @assert L > 0 "Vehicle length ($L) must be greater than 0"
        @assert L > lr "Vehicle length ($L) must be greater than the distance to the rear wheel ($lr)"
        new(ref, L, lr)
    end
end

function dynamics(model::BicycleModel, x, u)
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
    return SA[ẋ, ẏ, ω, ϕ]
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
