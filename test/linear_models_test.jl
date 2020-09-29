using RobotZoo
using RobotDynamics
using StaticArrays
using BenchmarkTools

linmodel = RobotZoo.LinearModels.DoubleIntegrator(3)
A = RobotDynamics.get_A(linmodel)
B = RobotDynamics.get_B(linmodel)
@test A isa SizedMatrix{6,6}
@test B isa SizedMatrix{6,3}
@test !RobotDynamics.is_discrete(linmodel)

model = RobotZoo.DoubleIntegrator(3)

F1 = RobotDynamics.DynamicsJacobian(linmodel)
F2 = RobotDynamics.DynamicsJacobian(model)
z = KnotPoint(rand(model)..., 0.1)
b1 = @benchmark dynamics($linmodel, $z)
b2 = @benchmark dynamics($model, $z)
@test minimum(b2) < minimum(b1)  # standard model is faster since it avoids matrix mult

b1 = @benchmark jacobian!($F1, $linmodel, $z)
b2 = @benchmark jacobian!($F2, $model, $z)
@test minimum(b1) < minimum(b2)  # linear model Jacobian is faster