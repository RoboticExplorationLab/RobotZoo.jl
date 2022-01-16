module LinearModels

using RobotDynamics

"""
    DoubleIntegrator(D=1)

Create a LinearModel of a `D`-dimensional double integrator, such that the state dimension 
is `2*D` (`D` positions and `D` velocities). The model is continuous and non-affine.

This is, in practice, identical to the `DoubleIntegrator` type in the base `RobotZoo`; 
however the other is a child of `ContinuousDynamics` whereas this one is explicitly a 
`LinearModel`. Both are included for reference and performance comparisons.
""" 
function DoubleIntegrator(D=1)
    A = zeros(2D,2D)
    B = zeros(2D,D)
    for i = 1:D
        A[i, D+i] = 1
        B[D+i, i] = 1
    end
    return LinearModel(A, B)
end

end
