struct PlanarQuadrotor <: AbstractModel
    mass::Float64  # mass (kg)
    g::Float64     # gravity (m/s²)
    ℓ::Float64     # tip to tip distance (m)
    J::Float64
    function PlanarQuadrotor(mass::Real, g::Real, ℓ::Real)
        J = 0.2*mass*ℓ^2
        new(mass,g,ℓ,J)
    end
end
PlanarQuadrotor(; mass=1.0, g=9.81, ℓ=0.3) = PlanarQuadrotor(mass, g, ℓ)

state_dim(::PlanarQuadrotor) = 6
control_dim(::PlanarQuadrotor) = 2

function dynamics(model::PlanarQuadrotor, x, u)
    mass,g,ℓ,J = model.mass, model.g, model.ℓ, model.J

    s,c = sincos(θ)
    θ = x[3]
    ẍ = (1/mass)*(u[1] + u[2])*s
    ÿ = (1/mass)*(u[1] + u[2])*c - g
    θddot = (1/J)*(ℓ/2)*(u[2] - u[1])

    return SA[x[4], x[5], x[6], ẍ, ÿ, θddot]
end
