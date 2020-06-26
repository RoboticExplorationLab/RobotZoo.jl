@with_kw struct FreeBody{R,T} <: RigidBody{R}
    mass::T = 1.0
    J::Diagonal{T,SVector{3,T}} = Diagonal(@SVector ones(3))
    Jinv::Diagonal{T,SVector{3,T}} = Diagonal(@SVector ones(3))
end
RobotDynamics.control_dim(::FreeBody) = 6

(::FreeBody)(;kwargs...) = FreeBody{UnitQuaternion{Float64,CayleyMap},Float64}(;kwargs...)

function forces(model::FreeBody, x::StaticVector, u::StaticVector)
    q = orientation(model, x)
    F = @SVector [u[1], u[2], u[3]]
    q*F  # world frame
end

function moments(model::FreeBody, x::StaticVector, u::StaticVector)
    return @SVector [u[4], u[5], u[6]]  # body frame
end

function RobotDynamics.wrench_jacobian!(F, model::FreeBody, z)
    x = state(z)
    u = control(z)
    q = orientation(model, x)
    ir, iq, iv, iω, iu = RobotDynamics.gen_inds(model)
    iF = SA[1,2,3]
    iM = SA[4,5,6]
    F[iF, iq] .= Rotations.∇rotate(q, u[iF])
    F[iF, iu[iF]] .= RotMatrix(q)
    for i = 1:3
        F[iM[i], iu[i+3]] = 1
    end
    return F
end

function RobotDynamics.wrench_sparsity(::FreeBody)
    SA[false true  false false true;
       false false false false true]
end

inertia(model::FreeBody, x, u) = model.J
inertia_inv(model::FreeBody, x, u) = model.Jinv
mass(model::FreeBody, x, u) = model.mass
