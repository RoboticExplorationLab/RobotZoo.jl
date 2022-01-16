@autodiff struct PlanarRocket{T} <: ContinuousDynamics
    g::T  # gravity (m/s²)
    m::T  # mass (kg)
    ℓ::T  # length (m)
    max_angle::T  # max thrust angle (deg)
    max_roll::T   # max roll angle (deg)
    Tlo::T        # min thrust, per unit weight (1/kg⋅m/s²)
    Thi::T        # max thrust, per unit weight (1/kg⋅m/s²)
    function PlanarRocket(;g=9.81, m=100., ℓ=2.0, max_angle=10., max_roll=10., Tlo=0.75, Thi=1.5)
        # g,m,ℓ,max_angle,max_roll,Tlo,Thi = promote(g,m,ℓ,max_angle,max_roll,Tlo,Thi)
        T = eltype(promote(g,m,ℓ,max_angle,max_roll,Tlo,Thi))
        new{T}(g,m,ℓ,max_angle,max_roll,Tlo,Thi)
    end
end
state_dim(::PlanarRocket) = 8
control_dim(::PlanarRocket) = 2

umin(model::PlanarRocket) = SA[model.Tlo*model.m*model.g, -deg2rad(model.max_angle)]
umax(model::PlanarRocket) = SA[model.Thi*model.m*model.g, +deg2rad(model.max_angle)]

begin
    body = quote
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
    end
    @eval function dynamics(model::PlanarRocket, x, u)
        $body
        return SA[x[4], x[5], x[6], a[1], a[2], a[3], u[1], u[2]]
    end
    @eval function dynamics!(model::PlanarRocket, xdot, x, u)
        $body
        xdot[1:3] .= @view x[4:6]
        xdot[4:6] .= a
        xdot[7:8] .= u
        return nothing
    end
end