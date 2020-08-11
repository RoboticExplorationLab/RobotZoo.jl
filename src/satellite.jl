import Rotations: lmult, vmat, hmat

"""
    Satellite

A simple satellite model, modeling only the attitude dynamics with no external forces. Has
direct control over torques about the body axes.

# Constructor
    Satellite()
    Satellite(J)

where `J` is a `Diagonal` inertia matrix.
"""
struct Satellite <: LieGroupModel
    J::Diagonal{Float64,SVector{3,Float64}}
end

Satellite() = Satellite(Diagonal(@SVector ones(3)))

RobotDynamics.control_dim(::Satellite) = 3
Base.position(::Satellite, x) = @SVector zeros(3)
orientation(::Satellite, x) = UnitQuaternion(x[4], x[5], x[6], x[7])

RobotDynamics.LieState(::Satellite) = RobotDynamics.LieState(UnitQuaternion{Float64}, (3,0))

function dynamics(model::Satellite, x::SVector, u::SVector)
    ω = @SVector [x[1], x[2], x[3]]
    q = normalize(@SVector [x[4], x[5], x[6], x[7]])
    J = model.J

    ωdot = J\(u - ω × (J*ω))
    qdot = 0.5*lmult(q)*hmat()*ω
    return [ωdot; qdot]
end
