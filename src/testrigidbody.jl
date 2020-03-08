abstract type TestRigidBody{B,R} <: RigidBody{R} end

############################################################################################
#                            OVERRIDE ALL ROTATION STUFF
############################################################################################

Dynamics.state_diff_size(model::TestRigidBody{false,<:UnitQuaternion}) = 13
Dynamics.state_diff_size(model::TestRigidBody{false}) = 12

function Dynamics.state_diff(model::TestRigidBody{false}, x::SVector, x0::SVector)
    x - x0
end

function Dynamics.state_diff_jacobian(model::TestRigidBody{false}, x::SVector)
    return I
end

function ∇²differential(model::TestRigidBody{false},
        x::SVector, dx::SVector)
    return I*0
end

function inverse_map_jacobian(model::TestRigidBody{false}, x::SVector)
    return I
end

function TrajectoryOptimization.inverse_map_∇jacobian(model::TestRigidBody{false},
        x::SVector, b::SVector)
    return I*0
end
