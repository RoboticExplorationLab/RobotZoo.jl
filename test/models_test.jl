function num_allocs(model)
    dt = 0.1
    x, u = rand(model)
    n,m = size(model)
    z = KnotPoint(x, u, dt)
    ∇c  = zeros(n,n+m)
    dynamics(model, x, u)
    jacobian!(∇c, model, z)
    discrete_dynamics(RK3, model, x, u, z.t, dt)
    discrete_jacobian!(RK3, ∇c, model, z)
    allocs  = @allocated dynamics(model, x, u)
    allocs += @allocated jacobian!(∇c, model, z)
    allocs += @allocated discrete_dynamics(RK3, model, x, u, z.t, dt)
    allocs += @allocated discrete_jacobian!(RK3, ∇c, model, z)
end

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

# Car
car = RobotZoo.DubinsCar()
@test size(car) == (3,2)
@test num_allocs(car) == 0

# Cartpole
cartpole = RobotZoo.Cartpole()
@test size(cartpole) == (4,1)
@test_broken num_allocs(cartpole) == 0

# Quadrotor
quad = RobotZoo.Quadrotor()
@test size(quad) == (13,4)
@test_broken num_allocs(cartpole) == 0


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
