const RD = RobotDynamics
using RobotDynamics: KnotPoint, dynamics, dynamics!, jacobian!
using RobotDynamics: StaticReturn, InPlace, ForwardAD, FiniteDifference
using Random
function test_model(model; evals=1, samples=1, tol=1e-6) 
    println(typeof(model))
    dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
    t,dt = 1.1,0.1
    Random.seed!(1)
    x, u = rand(model)
    n,m = RD.dims(model)
    z = KnotPoint(x, u, t, dt)
    ∇c1  = zeros(n,n+m)
    ∇c2  = zeros(n,n+m)
    xdot = zeros(n)
    allocs = 0
    allocs += @ballocated RobotDynamics.dynamics($model, $x, $u) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.dynamics!($model, $xdot, $x, $u) evals=evals samples=samples
    @test xdot == RobotDynamics.dynamics(model, x, u)
    @test allocs == 0 
    RobotDynamics.jacobian!(StaticReturn(), ForwardAD(), model, ∇c1, xdot, z)
    RobotDynamics.jacobian!(StaticReturn(), FiniteDifference(), model, ∇c2, xdot, z)
    @test ∇c1 ≈ ∇c2 atol=tol
    RobotDynamics.jacobian!(InPlace(), ForwardAD(), model, ∇c2, xdot, z)
    @test ∇c1 ≈ ∇c2
    RobotDynamics.jacobian!(InPlace(), FiniteDifference(), model, ∇c1, xdot, z)
    @test ∇c1 ≈ ∇c2 atol=tol
end

# Acrobot
acrobot = RobotZoo.Acrobot()
@test RD.dims(acrobot) == (4,1,4)
test_model(acrobot)

# Car
car = RobotZoo.DubinsCar()
@test RD.dims(car) == (3,2,3)
test_model(car)

# Bicycle Car
bicycle = RobotZoo.BicycleModel()
@test RD.dims(bicycle) == (4,2,4)
test_model(bicycle)

bicycle = RobotZoo.BicycleModel(ref=:rear)
@test RD.dims(bicycle) == (4,2,4)
test_model(bicycle)

# Planar Rocket
rocket = RobotZoo.PlanarRocket()
@test RD.dims(rocket) == (8,2,8)
test_model(rocket)

# Planar Quad
quad = RobotZoo.PlanarQuadrotor()
@test RD.dims(quad) == (6,2,6)
test_model(quad)

# Cartpole
cartpole = RobotZoo.Cartpole()
@test RD.dims(cartpole) == (4,1,4)
test_model(cartpole)

# Double Integrator
dim = 3
di = RobotZoo.DoubleIntegrator(dim)
n,m = RD.dims(di)
@test (n,m) == (6,3)
test_model(di)

# Pendulum
pend = RobotZoo.Pendulum()
RobotZoo.Pendulum{Float64}(1,1,1, 1,1,1)
@test RD.dims(pend) == (2,1,2)
test_model(pend)

# Quadrotor
quad = RobotZoo.Quadrotor()
@test RD.dims(quad) == (13,4,13)
test_model(quad,tol=1e-4)

# Yak Plane
yak = RobotZoo.YakPlane(MRP{Float64})
@test RD.dims(yak) == (12,4,12)
test_model(yak,tol=1e-4)

# Test other functions
dt = 0.1
car = RobotZoo.DubinsCar()
n,m = RD.dims(car)
@test zeros(car) == (zeros(n), zeros(m))
@test zeros(Int,car)[1] isa SVector{n,Int}
@test fill(car,0.1) == (fill(0.1,n), fill(0.1,m))
# @test ones(Float32,car)[2] isa SVector{m,Float32}
