using RobotZoo
using RobotDynamics
using StaticArrays
using BenchmarkTools
using FiniteDiff
using LinearAlgebra
using Test
using Random
Random.seed!(1)
const RD = RobotDynamics

# Continuous Linear model
n,m = 10,5
A,B = RobotZoo.RandomLinearModels.gencontrollable(n,m, :continuous)
d = randn(10)
model = RobotZoo.LinearModel(A,B,d)
@test RD.dims(model) == (n,m,n)
x,u = rand(model)
@test length(x) == n
@test length(u) == m
@test RD.dynamics(model, x, u) ≈ A*x + B*u + d
xn = similar(x)
RD.dynamics!(model, xn, x, u)
@test xn ≈ A*x + B*u + d

@test RobotZoo.iscontrollable(model)

J = zeros(n, n+m)
z = RD.KnotPoint(x, u, 0, NaN)
RD.jacobian!(RD.InPlace(), RD.UserDefined(), model, J, xn, z)
@test J[:,1:n] ≈ A
@test J[:,n+1:end] ≈ B

Jfd = similar(J)
FiniteDiff.finite_difference_jacobian!(
    Jfd, 
    (y,z)->RD.dynamics!(model, y, z[1:n], z[n+1:end]), 
    Vector(z.z)
)
@test Jfd ≈ J rtol=1e-6

# Linearized Continuous Model
Random.seed!(1)
model0 = RobotZoo.Cartpole()
xe = [0,pi,0,0]
ue = [0.]
model = RobotZoo.LinearModel(RD.InPlace(), RD.ForwardAD(), model0, xe, ue)
n,m = RD.dims(model)
@test n == 4
@test m == 1
dx = randn(4)
du = randn(1)
@test RD.dynamics(model, dx*0, du*0) ≈ RD.dynamics(model0, xe, ue)
@test norm(RD.dynamics(model, dx*1e-3, du*1e-3) - RD.dynamics(model0, xe, ue)) < 0.1

Jfd = zeros(4,5) 
y = zeros(4)
z = RD.KnotPoint{n,m}(xe, ue, 0, NaN)
FiniteDiff.finite_difference_jacobian!(
    Jfd, 
    (y,z)->RD.dynamics!(model0, y, z[1:n], z[n+1:end]), 
    Vector(z.z)
)
@test Jfd[:,1:n] ≈ RobotZoo.getstatecoeffs(model)
@test Jfd[:,n+1:end] ≈ RobotZoo.getcontrolcoeffs(model)
@test RobotZoo.isstable(model) == false

# Linearize about stable equilibrium
xe = [0,0.0,0,0]
ue = [0.]
model_stable = RobotZoo.LinearModel(RD.InPlace(), RD.ForwardAD(), model0, xe, ue)
@test RobotZoo.isstable(model_stable)

# Use random generator on linear model
n,m = 12,6
model_rand = rand(RobotZoo.LinearModel, n, m)
@test RD.dims(model_rand) == (n,m,n)
@test RobotZoo.iscontrollable(model_rand)

## Discrete Linear Model
n,m = 10,5
h = 0.1
A,B = RobotZoo.RandomLinearModels.gencontrollable(n, m, :discrete)
d = randn(n)
model = RobotZoo.DiscreteLinearModel(h,A,B)
@test RobotZoo.getstatecoeffs(model) == A 
@test RobotZoo.getcontrolcoeffs(model) == B 
@test RobotZoo.getaffinecoeffs(model) == zeros(n)
@test RobotZoo.iscontrollable(model)

model = RobotZoo.DiscreteLinearModel(h,A,B,d)
@test RobotZoo.getstatecoeffs(model) == A
@test RobotZoo.getcontrolcoeffs(model) == B
@test RobotZoo.getaffinecoeffs(model) == d 
@test RobotZoo.iscontrollable(model)

C = randn(n,n)
model = RobotZoo.DiscreteLinearModel(h,A,B,C,d)
@test RobotZoo.getstatecoeffs(model) == A
@test RobotZoo.getcontrolcoeffs(model) == B
@test RobotZoo.getnextstatecoeffs(model) == C
@test RobotZoo.getaffinecoeffs(model) == d 

@test_nowarn RobotZoo.checktimestep(model, h)
@test_nowarn RobotZoo.checktimestep(model, h + eps(h))
@test_throws ArgumentError RobotZoo.checktimestep(model, h*1.1)

@test RD.dims(model) == (n,m,n)
x,u = rand(model)
z = RD.KnotPoint(x,u,0.,h)

# Explicit
model_explicit = RobotZoo.DiscreteLinearModel(h,A,B,d)
@test RobotZoo.getnextstatecoeffs(model_explicit) ≈ -I(n)
@test RD.discrete_dynamics(model_explicit, x, u, 0, h) ≈ A*x + B*u + d 

J = zeros(n,n+m)
xn = zeros(n)
RD.jacobian!(RD.InPlace(), RD.UserDefined(), model_explicit, J, xn, z) 
@test J[:,1:n] ≈ A
@test J[:,n+1:end] ≈ B

# Implicit
model_implicit = RobotZoo.DiscreteLinearModel(h,A,B,C,d)
@test RobotZoo.getnextstatecoeffs(model_implicit) ≈ C 
@test RD.discrete_dynamics(model_implicit, x, u, 0, h) ≈ -C\(A*x + B*u + d)

RD.jacobian!(RD.InPlace(), RD.UserDefined(), model_implicit, J, xn, z) 
@test J[:,1:n] ≈ -C\A
@test J[:,n+1:end] ≈ -C\B

## Discretize and Linearize
model0 = RobotZoo.Cartpole()
n,m = RD.dims(model0)

# Implicit
dmodel_explicit = RD.DiscretizedDynamics{RD.RK3}(model0)
xe = [0,pi,0,0]
ue = [0.]
model_explicit = RobotZoo.DiscreteLinearModel(
    RD.InPlace(),
    RD.ForwardAD(),
    dmodel_explicit,
    xe, ue, 0., h
)

J = zeros(n, n+m)
z = RD.KnotPoint{n,m}([xe;ue], 0.0, h)
xn = zeros(n)
RD.jacobian!(RD.StaticReturn(), RD.ForwardAD(), dmodel_explicit, J, xn, z)
@test J[:,1:n] ≈ RobotZoo.getstatecoeffs(model_explicit)
@test J[:,n+1:end] ≈ RobotZoo.getcontrolcoeffs(model_explicit)
@test RD.discrete_dynamics(dmodel_explicit, z) ≈ 
    RobotZoo.getaffinecoeffs(model_explicit)
@test RobotZoo.getnextstatecoeffs(model_explicit) == -I(n) * one(eltype(A))
@test RobotZoo.isstable(model_explicit) == false
@test RobotZoo.iscontrollable(model_explicit)

xe_stable = zeros(n)
model_explicit_stable = RobotZoo.DiscreteLinearModel(
    RD.InPlace(),
    RD.ForwardAD(),
    dmodel_explicit,
    xe_stable, ue, 0., h
)
@test RobotZoo.isstable(model_explicit_stable) == true 

# Implicit
dmodel_implicit = RD.DiscretizedDynamics{RD.ImplicitMidpoint}(model0)
model_implicit= RobotZoo.DiscreteLinearModel(
    RD.InPlace(),
    RD.ForwardAD(),
    dmodel_implicit,
    xe, ue, 0., h
)
@test RD.dims(model_implicit) == (n,m,n)

Jc = zeros(n,n+m)
RD.jacobian!(RD.StaticReturn(), RD.ForwardAD(), model0, Jc, xn, z)
A = Jc[:,1:n]
B = Jc[:,n+1:end]
@test RobotZoo.getstatecoeffs(model_implicit) ≈ h/2*A + I
@test RobotZoo.getnextstatecoeffs(model_implicit) ≈ h/2*A - I
@test RobotZoo.getcontrolcoeffs(model_implicit) ≈ h*B

# Exponential integrator
model_linear = RobotZoo.LinearModel(RD.InPlace(), RD.ForwardAD(), 
    model0, xe, ue, affine=false)
model_exp = RobotZoo.DiscreteLinearModel(model_linear, h)

# Compare the linear models
dx = randn(n) * 1e-2
du = randn(m) * 1e-2
xn_exp = RD.discrete_dynamics(model_exp, dx, du, 0.0, h)
xn_implicit = RD.discrete_dynamics(model_implicit, dx, du, 0, h) - model_implicit.d
xn_explicit = RD.discrete_dynamics(model_explicit, dx, du, 0, h) - model_explicit.d
@test norm(xn_exp - xn_implicit) < 1e-3
@test norm(xn_exp - xn_explicit) < 1e-3
