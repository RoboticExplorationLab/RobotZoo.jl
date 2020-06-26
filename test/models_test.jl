function num_allocs(model; evals=10, samples=2)
    dt = 0.1
    x, u = rand(model)
    n,m = size(model)
    z = KnotPoint(x, u, dt)
    ∇c  = zeros(n,n+m)
    b1 = @benchmark dynamics($model, $x, $u) evals=evals samples=samples
    b2 = @benchmark jacobian!($∇c, $model, $z) evals=evals samples=samples
    b3 = @benchmark discrete_dynamics($RK3, $model, $x, $u, $z.t, $dt) evals=evals samples=samples
    b4 = @benchmark discrete_jacobian!($RK3, $∇c, $model, $z) evals=evals samples=samples
    return b1.allocs + b2.allocs + b3.allocs + b4.allocs
end

# Acrobot
acrobot = RobotZoo.Acrobot()
@test size(acrobot) == (4,1)
@test num_allocs(acrobot) == 0

# Car
car = RobotZoo.DubinsCar()
@test size(car) == (3,2)
@test num_allocs(car) == 0

# Cartpole
cartpole = RobotZoo.Cartpole()
@test size(cartpole) == (4,1)
@test num_allocs(cartpole) == 0

# Double Integrator
dim = 3
di = RobotZoo.DoubleIntegrator(dim)
n,m = size(di)
@test (n,m) == (6,3)
x,u = rand(di)
@test num_allocs(di) == 0

# Pendulum
pend = RobotZoo.Pendulum()
@test size(pend) == (2,1)
@test num_allocs(pend) == 0

# Quadrotor
quad = RobotZoo.Quadrotor()
@test size(quad) == (13,4)
@test num_allocs(quad) == 0

# Yak Plane
yak = RobotZoo.YakPlane(MRP{Float64})
@test size(yak) == (12,4)
@test num_allocs(yak) == 0

# Test other functions
dt = 0.1
car = RobotZoo.DubinsCar()
n,m = size(car)
@test zeros(car) == (zeros(n), zeros(m))
@test zeros(Int,car)[1] isa SVector{n,Int}
@test fill(car,0.1) == (fill(0.1,n), fill(0.1,m))
@test ones(Float32,car)[2] isa SVector{m,Float32}

# Test default integrator
x,u = rand(car)
z = KnotPoint(x,u,dt)
@test discrete_dynamics(car, z) == discrete_dynamics(RK3, car, z)
