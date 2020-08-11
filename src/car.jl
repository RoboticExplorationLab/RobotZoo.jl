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
