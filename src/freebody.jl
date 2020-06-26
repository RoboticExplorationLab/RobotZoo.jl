"""
    FreeBody{R,T}

A rigid body moving freely in 3-dimensional space, with direct control of both force and torque
at the center of mass, for a total of 6 controls.

# Constructors
    FreeBody{R,T}(mass, J)
    FreeBody(R=UnitQuaternion{Float64}; kwargs...)

where `R <: Rotation{3}` with keyword arguments
* `mass` - mass of the body
* `J` - a vector of the three diagonal elements of the inertia matrix
"""
struct FreeBody{R,T} <: RigidBody{R}
    mass::T
    J::Diagonal{T,SVector{3,T}}
    Jinv::Diagonal{T,SVector{3,T}}
    function FreeBody{R,T}(mass::Real, Jdiag::AbstractVector) where {R,T}
        @assert length(Jdiag) == 3 "inertia diagonal must have only 3 entries"
        @assert mass > 0 "mass must be positive"
        J = SVector{3}(Jdiag)
        new{R,T}(mass, J, inv(J))
    end
end
RobotDynamics.control_dim(::FreeBody) = 6

FreeBody(R=UnitQuaternion{Float64}; mass=1.0, J=@SVector ones(3)) =
    FreeBody{R,promote_type(typeof(mass), eltype(J))}(mass, J)

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
