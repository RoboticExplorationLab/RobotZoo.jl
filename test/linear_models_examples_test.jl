using RobotZoo
using RobotDynamics
using StaticArrays
using BenchmarkTools
using FiniteDiff
using LinearAlgebra
using Test
using Random

Random.seed!(2)
const RD = RobotDynamics

## Double Integrator
D = 2
linmodel = RobotZoo.LinearModels.DoubleIntegrator(D)
model = RobotZoo.DoubleIntegrator(D)
x, u = rand(model)
@test RD.dynamics(model, x, u) ≈ RD.dynamics(linmodel, x, u)

h = 0.1
dmodel_exp = RobotZoo.DiscreteLinearModel(linmodel, h)
dmodel_lin = RobotZoo.LinearModels.DiscreteDoubleIntegrator(h, D)
@test RD.discrete_dynamics(dmodel_exp, x, u, 0.0, h) ≈
      RD.discrete_dynamics(dmodel_lin, x, u, 0.0, h)

## FlexibleSatellite
model = RobotZoo.LinearModels.FlexibleSatellite()
@test RD.dims(model) == (12, 3, 12)