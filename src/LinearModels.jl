module LinearModels

using RobotDynamics
using LinearAlgebra
import RobotZoo: LinearModel, DiscreteLinearModel

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

function DiscreteDoubleIntegrator(h,D=1)
    b = h^2 / 2
    A = zeros(2D, 2D)
    B = zeros(2D, D)
    for i = 1:2D
        A[i,i] = 1.0
    end
    for i = 1:D
        A[i, i+D] = h
        B[i, i] = b
        B[i+D, i] = h
    end
    DiscreteLinearModel(h, A, B)
end

function FlexibleSatellite()
    eye(n) = Matrix(I,n,n) 

    # inertia matrix
    J = diagm([1;2;3])

    # reaction wheel jacobian
    B_sc = diagm(ones(3))


    # linear momentum coupling matrix
    phi = [0 1 0;
           1 0 0;
           0 .2 -.8];

    # angular momentum coupling matrix
    delta = [0 0 1;
             0 1 0;
            -.7 .1 .1]

    # store this matrix for faster computations
    T = inv(J-delta'*delta)

    j = 3; # 3 modes

    # damping and stiffness
    zeta = [.001;.001;.001]
    Delta = [.05; .2; .125] * (2*pi)

    # damping and stiffness matrices
    C = zeros(j,j)
    K = zeros(j,j)
    for i =1:j
        C[i,i] = 2*zeta[i]*Delta[i];
        K[i,i] = Delta[i]^2;
    end


               #   mrp        w                  n                       ndot
    pdot_row = [zeros(3,3) .25*eye(3)       zeros(3,j)                 zeros(3,j)];
    wdot_row = [zeros(3,3) zeros(3,3)     T*delta'*K                  T*delta'*C];
    ndot_row = [zeros(j,3) zeros(j,3)     zeros(j,j)                  eye(j)];
    nddot_row = [zeros(j,3) zeros(j,3) (-K - delta*T*delta'*K)    (-C - delta*T*delta'*C)];

    # analytical A
    A_analytical = [pdot_row;wdot_row;ndot_row;nddot_row];

    # analytical B
    B_analytical = [zeros(3,3);
              -T*B_sc;
              zeros(j,3);
              delta*T*B_sc];

    LinearModel(A_analytical, B_analytical)
end

end