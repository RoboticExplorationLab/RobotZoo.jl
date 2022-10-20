
"""
    Crazyflie{R}

A crazyflie model, with simple aerodynamic forces. The orientation is represent by
a general rotation `R`. The body z-axis point is vertical, so positive controls cause acceleration
in the positive z direction.

# Constructor
    Crazyflie(; kwargs...)
    Crazyflie{R}(; kwargs...)

where `R <: Rotation{3}` and defaults to `UnitQuaternion{Float64}` if omitted. The keyword arguments are
* `mass` - mass of the crazyflie, in kg (default = 0.027)
* `J` - inertia of the crazyflie, in kg⋅m² (default = `Diagonal([0.00003144988, 0.00003151127, 0.00007058874])`)
* `gravity` - gravity vector, in kg/m² (default = [0,0,-9.81])
* `motor_dist` - distane between the motors, in m (default = 0.065)
* `km` - motor torque constant (default = 0.0000024)
* `kf` - motor force constant (default = (0.06*9.81/4) / 65536)
"""
@autodiff struct Crazyflie{R} <: RigidBody{R}
    mass::Float64
    J::SMatrix{3,3,Float64,9}
    Jinv::SMatrix{3,3,Float64,9}
    gravity::SVector{3,Float64}
    motor_dist::Float64
    kf::Float64
    km::Float64
    bodyframe::Bool  # velocity in body frame?
    ned::Bool
end
control_dim(::Crazyflie) = 4

function Crazyflie{R}(;
        mass = 0.5,
        J = Diagonal(@SVector [0.00003144988, 0.00003151127, 0.00007058874]),
        gravity = SVector(0, 0, -9.81),
        motor_dist = 0.065 / 2,
        kf = (0.06*9.81/4) / 65536,
        km = 0.005964552,
        bodyframe = false,
        ned = false,
    ) where R
    @assert issymmetric(J)
    Crazyflie{R}(mass, J, inv(J), gravity, motor_dist, kf, km, bodyframe, ned)
end

(::Type{Crazyflie})(;kwargs...) = Crazyflie{UnitQuaternion{Float64}}(;kwargs...)

@inline RobotDynamics.velocity_frame(model::Crazyflie) = model.bodyframe ? :body : :world

function trim_controls(model::Crazyflie)
    @SVector fill(-model.gravity[3]*model.mass/4.0, size(model)[2])
end

function forces(model::Crazyflie, x, u)
    q = orientation(model, x)
    kf = model.kf
    g = model.gravity
    m = model.mass

    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]

    F1 = max(0,kf*w1);
    F2 = max(0,kf*w2);
    F3 = max(0,kf*w3);
    F4 = max(0,kf*w4);
    F = @SVector [0., 0., F1+F2+F3+F4] #total rotor force in body frame
    if model.ned
        F = SA[0,0,-F[3]]
        g = -g
    end

    f = m*g + q*F # forces in world frame
    return f
end

function moments(model::Crazyflie, x, u)

    kf, km = model.kf, model.km
    L = model.motor_dist

    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]

    F1 = max(0,kf*w1);
    F2 = max(0,kf*w2);
    F3 = max(0,kf*w3);
    F4 = max(0,kf*w4);

    M1 = km*w1;
    M2 = km*w2;
    M3 = km*w3;
    M4 = km*w4;
    tau = @SVector [L*(F3+F4-F1-F2), L*(F1+F4-F2-F3), (M1-M2+M3-M4)] #total rotor torque in body frame
    if model.ned
        tau = SA[tau[1], -tau[2], -tau[3]]
    end
    return tau
end

function wrenches(model::Crazyflie, x, u)
    F = forces(model, x, u)
    M = moments(model, x, u)
    return [F; M]

    q = orientation(model, x)
    C = forceMatrix(model)
    mass, g = model.mass, model.gravity

    # Calculate force and moments
    w = max.(u, 0)  # keep forces positive
    fM = forceMatrix(model)*w
    f = fM[1]
    M = @SVector [fM[2], fM[3], fM[4]]
    e3 = @SVector [0,0,1]
    F = mass*g - q*(f*e3)
    return F,M
end

function forceMatrix(model::Crazyflie)
    kf, km = model.kf, model.km
    L = model.motor_dist
    @SMatrix [
        kf      kf      kf      kf;
        -L*kf   -L*kf   L*kf    L*kf;
        L*kf    -L*kf   -L*kf   L*kf;
        km      -km     km      -km;
    ]
end


RobotDynamics.inertia(model::Crazyflie) = model.J
RobotDynamics.inertia_inv(model::Crazyflie) = model.Jinv
RobotDynamics.mass(model::Crazyflie) = model.mass

function Base.zeros(model::Crazyflie{R}) where R
    x = RobotDynamics.build_state(model, zero(RBState))
    u = @SVector fill(-model.mass*model.gravity[end]/4, 4)
    return x,u
end
