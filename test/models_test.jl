const RD = RobotDynamics
using RobotDynamics: KnotPoint, dynamics, dynamics!, jacobian!
using RobotDynamics: StaticReturn, InPlace, ForwardAD, FiniteDifference
function num_allocs(model; evals=1, samples=1, inplace=true) 
    dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
    t,dt = 1.1,0.1
    x, u = rand(model)
    n,m = RD.dims(model)
    xdot = zeros(n)
    z = KnotPoint(x, u, t, dt)
    ∇c  = zeros(n,n+m)
    xn = zeros(n)
    allocs = 0
    allocs += @ballocated RobotDynamics.dynamics($model, $x, $u) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.dynamics!($model, $xdot, $x, $u) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.jacobian!(StaticReturn(), ForwardAD(), $model, $∇c, $xdot, $z) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.jacobian!(StaticReturn(), FiniteDifference(), $model, $∇c, $xdot, $z) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.jacobian!(InPlace(), ForwardAD(), $model, $∇c, $xdot, $z) evals=evals samples=samples
    allocs += @ballocated RobotDynamics.jacobian!(InPlace(), FiniteDifference(), $model, $∇c, $xdot, $z) evals=evals samples=samples

    if inplace
        allocs += @ballocated RobotDynamics.discrete_dynamics($dmodel, $z) evals=evals samples=samples
        allocs += @ballocated RobotDynamics.discrete_dynamics!($dmodel, $xdot, $z) evals=evals samples=samples
        allocs += @ballocated RobotDynamics.jacobian!(StaticReturn(), ForwardAD(), $dmodel, $∇c, $xdot, $z) evals=evals samples=samples
        allocs += @ballocated RobotDynamics.jacobian!(StaticReturn(), FiniteDifference(), $dmodel, $∇c, $xdot, $z) evals=evals samples=samples
        allocs += @ballocated RobotDynamics.jacobian!(InPlace(), ForwardAD(), $dmodel, $∇c, $xdot, $z) evals=evals samples=samples
        allocs += @ballocated RobotDynamics.jacobian!(InPlace(), FiniteDifference(), $dmodel, $∇c, $xdot, $z) evals=evals samples=samples
    end

    return allocs
end

# Acrobot
acrobot = RobotZoo.Acrobot()
@test RD.dims(acrobot) == (4,1,4)
@test num_allocs(acrobot) == 0

# Car
car = RobotZoo.DubinsCar()
@test RD.dims(car) == (3,2,3)
@test num_allocs(car) == 0

# Bicycle Car
bicycle = RobotZoo.BicycleModel()
@test RD.dims(bicycle) == (4,2,4)
@test num_allocs(bicycle) == 0

bicycle = RobotZoo.BicycleModel(ref=:rear)
@test RD.dims(bicycle) == (4,2,4)
@test num_allocs(bicycle) == 0

# Planar Rocket
rocket = RobotZoo.PlanarRocket()
@test RD.dims(rocket) == (8,2,8)
@test num_allocs(rocket) == 0

# Planar Quad
quad = RobotZoo.PlanarQuadrotor()
@test RD.dims(quad) == (6,2,6)
@test num_allocs(rocket) == 0

# Cartpole
cartpole = RobotZoo.Cartpole()
@test RD.dims(cartpole) == (4,1,4)
@test num_allocs(cartpole) == 0

# Double Integrator
dim = 3
di = RobotZoo.DoubleIntegrator(dim)
n,m = size(di)
@test (n,m) == (6,3)
@test num_allocs(di) == 0

# Pendulum
pend = RobotZoo.Pendulum()
RobotZoo.Pendulum{Float64}(1,1,1, 1,1,1)
@test RD.dims(pend) == (2,1,2)
@test num_allocs(pend) == 0

# Quadrotor
quad = RobotZoo.Quadrotor()
@test RD.dims(quad) == (13,4,13)
@test num_allocs(quad, inplace=false) == 0

# Yak Plane
yak = RobotZoo.YakPlane(MRP{Float64})
@test RD.dims(yak) == (12,4,12)
@test num_allocs(yak, inplace=false) == 0

# Test other functions
dt = 0.1
car = RobotZoo.DubinsCar()
n,m = size(car)
@test zeros(car) == (zeros(n), zeros(m))
@test zeros(Int,car)[1] isa SVector{n,Int}
@test fill(car,0.1) == (fill(0.1,n), fill(0.1,m))
# @test ones(Float32,car)[2] isa SVector{m,Float32}
