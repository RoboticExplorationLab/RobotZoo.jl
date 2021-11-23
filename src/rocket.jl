Base.@kwdef struct PlanarRocket{T} <: ContinuousDynamics
    g::T = 9.81   # gravity (m/s²)
    m::T = 100.0  # mass (kg)
    ℓ::T = 2.0    # length (m)
    max_angle::T = 10.0  # max thrust angle (deg)
    max_roll::T = 10.0   # max roll angle (deg)
    Tlo::T = 0.75        # min thrust, per unit weight (1/kg⋅m/s²)
    Thi::T = 1.50        # max thrust, per unit weight (1/kg⋅m/s²)
end
state_dim(::PlanarRocket) = 8
control_dim(::PlanarRocket) = 2

umin(model::PlanarRocket) = SA[model.Tlo*model.m*model.g, -deg2rad(model.max_angle)]
umax(model::PlanarRocket) = SA[model.Thi*model.m*model.g, +deg2rad(model.max_angle)]

function dynamics(model::PlanarRocket, x, u)
    mass,g,ℓ = model.m, model.g, model.ℓ

    # clamp the controls
    # u = clamp.(u, umin(model), umax(model))

    ulo = umin(model)
    uhi = umax(model)
    T = clamp(x[7],ulo[1],uhi[1])   # thrust (N)
    ϕ = clamp(x[8],ulo[2],uhi[2])  # gimbal angle (rad)

    J = mass*ℓ^2 / 12  # moment of inertia 
    θ = x[3]

    # acceleration
    s,c = sincos(θ + ϕ)
    a = SA[
        (T/mass) * s,
        (T/mass) * c - g,
        (T/J) * 0.5 * ℓ * sin(ϕ)
    ]
    return SA[x[4], x[5], x[6], a[1], a[2], a[3], u[1], u[2]]
end