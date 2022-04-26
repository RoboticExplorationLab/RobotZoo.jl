"""
    LinearModel{M,V}
    
A continuous-time dynamics model of the form:

```math
\\dot{x} = A x + B u + d
```
"""
struct LinearModel{M,V} <: RD.ContinuousDynamics
    A::M
    B::M
    d::V
    function LinearModel(A::M, B::M, d::V) where {M<:AbstractMatrix, V<:AbstractVector}
        if !(size(A,1) == size(B,1) == size(d,1))
            throw(DimensionMismatch("All inputs must have the same number of rows."))
        elseif size(d,2) != 1
            throw(DimensionMismatch("d must have a single column."))
        end
        new{M,V}(A,B,d)
    end
end

function LinearModel(A::AbstractMatrix, B::AbstractMatrix)
    d = similar(A, size(A,1))
    d .= 0
    LinearModel(A, B, d)
end

function LinearModel(sig::RD.FunctionSignature, diffmethod::RD.DiffMethod, 
                     model::RD.ContinuousDynamics, x, u, t=zero(eltype(x)); affine::Bool=true)
    n,m = RD.dims(model)
    J = zeros(n,n+m)
    d = zeros(n)
    z = RD.KnotPoint{n,m}(x, u, t, NaN)
    RD.jacobian!(sig, diffmethod, model, J, d, z) 
    if affine
        if sig == RD.InPlace()
            RD.dynamics!(model, d, x, u, t)
        else
            d .= RD.dynamics(model, x, u, t)
        end
    else
        d .= 0
    end
    A = J[:,1:n]
    B = J[:,n+1:n+m]
    LinearModel(A, B, d)
end

# Generates a random controllable linear system
function Base.rand(::Type{LinearModel}, n::Integer, m::Integer)
    A,B = RandomLinearModels.gencontrollable(n, m, :continuous)
    d = randn(n)
    LinearModel(A, B, d)
end

RD.output_dim(model::LinearModel) = length(model.d)
RD.state_dim(model::LinearModel) = size(model.A, 2)
RD.control_dim(model::LinearModel) = size(model.B, 2)
RD.default_signature(model::LinearModel) = RD.InPlace()
RD.default_diffmethod(model::LinearModel) = RD.UserDefined()

getstatecoeffs(model::LinearModel) = model.A
getcontrolcoeffs(model::LinearModel) = model.B
getaffinecoeffs(model::LinearModel) = model.d

function iscontrollable(model::LinearModel)
    A = getstatecoeffs(model)
    B = getcontrolcoeffs(model)
    RandomLinearModels.iscontrollable(A, B)
end

isstable(model::LinearModel) = RandomLinearModels.isstablec(getstatecoeffs(model))

function dynamics(model::LinearModel, x, u)
    model.A * x .+ model.B * u .+ model.d
end

function dynamics!(model::LinearModel, xn, x, u)
    xn .= model.d
    mul!(xn, model.A, x, 1.0, 1.0)
    mul!(xn, model.B, u, 1.0, 1.0)
    nothing
end

function jacobian!(model::LinearModel, J, xn, x, u)
    n,m = RD.dims(model)
    J[:,1:n] .= model.A
    J[:,n+1:n+m] .= model.B
    nothing
end

"""
    DiscreteLinearModel{M,V}
    
A discrete-time dynamics model of the form:

```math
\\ A x_1 + B u_1 + C x_2 + d = 0
```
"""
struct DiscreteLinearModel{M,V} <: RD.DiscreteDynamics
    A::M
    B::M
    C::M
    d::V
    dt::Float64
    CA::M
    CB::M
    Cd::V
    function DiscreteLinearModel(dt, A::M, B::M, C::M, d::V) where {M<:AbstractMatrix, V<:AbstractVector}
        if !(size(A,1) == size(B,1) == size(C,1) == size(d,1))
            throw(DimensionMismatch("All inputs must have the same number of rows."))
        elseif size(A,2) != size(C,2)
            throw(DimensionMismatch("A and C must have the same number of columns."))
        elseif size(d,2) != 1
            throw(DimensionMismatch("d must have a single column."))
        end
        if C isa SMatrix
            F = lu(C)
        else
            F = factorize(C)
        end
        CA = -(F\A)
        CB = -(F\B)
        Cd = -(F\d)
        new{M,V}(A,B,C,d,dt, CA,CB,Cd)
    end
end

function DiscreteLinearModel(dt, A::AbstractMatrix, B::AbstractMatrix, 
                             d::AbstractVector=(similar(A, size(A,1)) .= zero(eltype(A))))
    C = similar(A) * zero(eltype(A))
    for i = 1:size(C,1)
        C[i,i] = -one(eltype(A))
    end
    DiscreteLinearModel(dt, A, B, C, d)
end

function DiscreteLinearModel(sig::RD.FunctionSignature, diffmethod::RD.DiffMethod, 
                             model::RD.DiscretizedDynamics, 
                             x, u, t, h; affine::Bool=true)
    n,m = RD.dims(model)
    J = zeros(n,n+m)
    d = zeros(n)
    z = RD.KnotPoint{n,m}([x;u], t, h)
    RD.jacobian!(sig, diffmethod, model, J, d, z) 
    if affine
        if sig == RD.InPlace()
            RD.discrete_dynamics!(model, d, x, u, t, h)
        else
            d .= RD.discrete_dynamics(model, x, u, t, h)
        end
    else
        d .= 0
    end
    A = J[:,1:n]
    B = J[:,n+1:n+m]
    DiscreteLinearModel(h, A, B, d)
end

function DiscreteLinearModel(sig::RD.FunctionSignature, diffmethod::RD.DiffMethod, 
                             model::RD.ImplicitDynamicsModel, 
                             x, u, t, h; affine::Bool=true)
    n,m = RD.dims(model)
    J1 = zeros(n,n+m)
    J2 = zeros(n,n+m)
    d = zeros(n)
    y1 = similar(d) 
    z1 = RD.StaticKnotPoint{Any,Any}(n,m, [x;u], t, h)
    z2 = RD.StaticKnotPoint{Any,Any}(n,m, [x;u], t, h)
    RD.dynamics_error_jacobian!(sig, diffmethod, model, J2, J1, d, y1, z2, z1) 
    if affine
        if sig == RD.InPlace()
            RD.dynamics_error!(model, d, y1, z2, z1)
        else
            d .= RD.dynamics_error(model, d, z2, z1)
        end
    else
        d .= 0
    end
    A = J1[:,1:n]
    B = J1[:,n+1:n+m]
    C = J2[:,1:n]
    if norm(J2[:,n+1:n+m]) > 0
        error("Cannot form a DiscreteLinearModel from an implicit dynamics function dependent on u2.")
    end
    DiscreteLinearModel(h, A, B, C, d)
end

# Exponential integration
function DiscreteLinearModel(model::LinearModel, h)
    if norm(getaffinecoeffs(model)) > 0
        throw(ArgumentError("Cannot use exponential integrator for an integrator with an affine term."))
    end
    n,m = RD.dims(model)
    A = getstatecoeffs(model)
    B = getcontrolcoeffs(model)
    expAB = exp([A*h B*h; zeros(m,n+m)])
    A = expAB[1:n,1:n]
    B = expAB[1:n,n+1:n+m]
    DiscreteLinearModel(h, A, B)
end

ouput_dim(model::DiscreteLinearModel) = size(model.A, 1)
state_dim(model::DiscreteLinearModel) = size(model.A, 2)
control_dim(model::DiscreteLinearModel) = size(model.B, 2)
RD.default_signature(model::DiscreteLinearModel) = RD.InPlace()
RD.default_diffmethod(model::DiscreteLinearModel) = RD.UserDefined()
getstatecoeffs(model::DiscreteLinearModel) = model.A
getcontrolcoeffs(model::DiscreteLinearModel) = model.B
getnextstatecoeffs(model::DiscreteLinearModel) = model.C
getaffinecoeffs(model::DiscreteLinearModel) = model.d

isstable(model::DiscreteLinearModel) = RandomLinearModels.isstabled(getstatecoeffs(model))

iscontrollable(model::DiscreteLinearModel) = RandomLinearModels.iscontrollable(model.CA, model.CB)

function checktimestep(model::DiscreteLinearModel, dt)
    if dt ≉ model.dt
        throw(ArgumentError("Invalid time step. Passed in $dt, but the model is defined using a time step of $(model.dt)."))
    end
end

function RD.discrete_dynamics(model::DiscreteLinearModel, x, u, t, h)
    checktimestep(model, h)
    model.CA * x + model.CB * u + model.Cd
end

function RD.discrete_dynamics!(model::DiscreteLinearModel, xn, x, u, t, h)
    checktimestep(model, h)
    xn .= model.Cd
    mul!(xn, model.CA, x, true, true)
    mul!(xn, model.CB, u, true, true) 
    nothing 
end

function RD.jacobian!(model::DiscreteLinearModel, J, xn, x, u, t, h)
    checktimestep(model, h) 
    n,m = RD.dims(model)
    J[:,1:n] .= model.CA
    J[:,n+1:n+m] .= model.CB
    nothing
end

function RD.dynamics_error(model::DiscreteLinearModel, 
                           z2::RD.AbstractKnotPoint, z1::RD.AbstractKnotPoint)
    checktimestep(model, RD.timestep(z1))
    x1,u1 = RD.state(z1), RD.control(z1)
    x2 = RD.state(z1)
    model.A * x1 + model.B * u1 + model.C * x2 + model.d
end

function RD.dynamics_error!(model::DiscreteLinearModel, y2, y1, 
                            z2::RD.AbstractKnotPoint, z1::RD.AbstractKnotPoint)
    checktimestep(model, RD.timestep(z1))
    x1,u1 = RD.state(z1), RD.control(z1)
    x2 = RD.state(z1)
    y2 .= model.d
    mul!(y2, model.A, x1, true, true)
    mul!(y2, model.B, u1, true, true)
    mul!(y2, model.C, x2, true, true)
    nothing
end

function RD.dynamics_error_jacobian!(model::DiscreteLinearModel, J2, J1, y2, y1, 
                                     z2::RD.AbstractKnotPoint, z1::RD.AbstractKnotPoint)
    checktimestep(model, RD.timestep(z1))
    n,m = RD.dims(model)
    J2[:,1:n] .= model.C
    J1[:,1:n] .= model.A
    J1[:,n+1:n+m] .= model.B
    nothing
end
