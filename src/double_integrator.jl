
struct DoubleIntegrator{N,M} <: AbstractModel
    pos::SVector{M,Int}
    vel::SVector{M,Int}
end

function DoubleIntegrator(D=1)
    pos = SVector{D,Int}(1:D)
    vel = SVector{D,Int}(D .+ (1:D))
    DoubleIntegrator{2D,D}(pos,vel)
end

RobotDynamics.state_dim(::DoubleIntegrator{N,M}) where {N,M} = N
RobotDynamics.control_dim(::DoubleIntegrator{N,M}) where {N,M} = M


@generated function dynamics(di::DoubleIntegrator{N,M}, x, u) where {N,M}
    vel = [:(x[$i]) for i = M+1:N]
    us = [:(u[$i]) for i = 1:M]
    :(SVector{$N}($(vel...),$(us...)))
end

Base.position(::DoubleIntegrator{<:Any,2}, x) = @SVector [x[1], x[2], 0]
orientation(::DoubleIntegrator, x) = UnitQuaternion(I)
