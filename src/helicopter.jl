import Rotations: lmult, vmat, hmat
using StaticArrays

@autodiff struct Helicopter{R} <: RigidBody{R}
    mass::Float64
    Ixx::Float64
    Iyy::Float64
    Izz::Float64
    J::SMatrix{3,3,Float64,9}
    Jinv::SMatrix{3,3,Float64,9}
    gravity::SVector{3,Float64}

    hmr::Float64
    omega_mr::Float64
    Rmr::Float64

    rho::Float64
    Cd::Float64
    Sx_fus::Float64
    Sy_fus::Float64
    Sz_fus::Float64
    
    ntr::Float64
    omega_tr::Float64
    Rtr::Float64
    ltr::Float64 
    htr::Float64 
 
end
control_dim(::Helicopter) = 4

function Helicopter{R}(;
    ########################################################################
    #### X-CELL 60 SE Helicopter Parameters: ###############################
    ########################################################################
    ########################################################################
    #=
    Vladislav Gavrilets. “Dynamic Model for a Miniature Aerobatic Helicopter”. 
    In: Handbook of Unmanned Aerial Vehicles. Ed. by Kimon P. Valavanis 
    and George J. Vachtsevanos. Dordrecht: Springer Netherlands, 2015, 
    pp. 279–306. ISBN : 978-90-481-9707-1. 
    DOI : 10 . 1007 / 978 - 90 - 481 - 9707 - 1 54. 
    URL : https : //doi.org/10.1007/978-90-481-9707-1 54.
    =#
    ########################################################################
    ## MASS
    mass = 8.2 #Kg

    ## INERTIA
    Ixx = 0.18 #kg m^2
    Iyy = 0.34 #kg m^2
    Izz = 0.28 #kg m^2
    J = Diagonal(@SVector[Ixx, Iyy, Izz])

    ## GRAVITY
    gravity = 9.81 #m/s^2

    ## MAIN ROTOR
    hmr = 0.235 #m
    omega_mr = 167 #rad/s
    Rmr = 0.775 #m

    ## AIR PARAMETERS
    rho = 1.204 #kg/m^3 #ρ

    ## DRAG PARAMETERS
    # https://www.engineeringtoolbox.com/drag-coefficient-d_627.html
    Cd = 1.17 #square flat plate 90degrees to the wind.
    Sx_fus = 0.1 #m^2
    Sy_fus = 0.22 #m^2
    Sz_fus = 0.15 #m^2
    
    ## TAIL ROTOR
    ntr = 4.66 #Gear ratio of tr to mr
    omega_tr = ntr * omega_mr #rad/s
    Rtr = 0.13 #m
    ltr = 0.91 #m
    htr = 0.08 #m 

    ) where R
    #@assert issymmetric(J)
    Helicopter{R}(mass,Ixx,Iyy,Izz,J,inv(J),hmr,omega_mr,Rmr,rho,Cd,Sx_fus,Sy_fus,Sz_fus,ntr,omega_tr,Rtr,ltr,htr)
end

#(::Type{Helicopter})(;kwargs...) = Helicopter{UnitQuaternion{Float64}}(;kwargs...)

#@inline RobotDynamics.velocity_frame(model::Quadrotor) = model.bodyframe ? :body : :world_vector


function trim_controls(model::Helicopter)
    @SVector [7.354261, 0.0685781, 80.4419975, -16.2049026] #Hover conditions found by root-finding algorithm
end

####################################################################
#### X-CELL 60 SE Helicopter Dynamics: #############################
####################################################################
####################################################################
function dynamics(model::Helicopter, x::SVector, u::SVector)
    """  
    #### NOTE ON STATES::
    x: [x y z Quat u v w p q r wx wy wz wxdot wydot wzdot]'                                                   
    x[14:16] --> wind velocities (estimated or ground-truth) 
    x[14:16] --> wind "accelerations" 

    #### NOTE ON CONTROL INPUTS:
    # u: [a b Tmr Ttr]'
    u[1] = Main Rotor tilt angle with vertical in Helicopter's X-Z plane (*check!) 
    u[2] = Main Rotor tilt angle with vertical in Helicopter's Y-Z plane
    u[3] = Main Rotor Thrust
    u[4] = Tail Rotor Thrust

    #### NOTE ON DYNAMICS:
    ASSUME: No blade-flapping dynamics;
            No horizontal+Vertical fins;
            No stabilizer;
            Centre of Pressure == C.O.G

    """

    m, g = model.mass, model.gravity
    hmr, Rmr, omega_mr = model.hmr, model.Rmr, model.omega_mr
    ltr, htr = model.ltr, model.htr
    rho = model.rho
    Cd  = model.Cd
    Sx_fus, Sy_fus, Sz_fus = model.Sx_fus, model.Sy_fus, model.Sz_fus

    r = SA[x[1], x[2], x[3]] #x y z
    q = SA[x[4], x[5], x[6], x[7]] # quaternion #NEED TO NORMAILIZE???
    v = SA[x[8], x[9], x[10]] #u v w
    ω = SA[x[11], x[12], x[13]] #p q r
    w = SA[x[14], x[15], x[16]] #wind-velocities
    vw = SA[x[17], x[18], x[19]] #wind-velocity_dot

    #Q = qtoQ(q) 
    Q = hmat()' * lmult(q) * rmult(q)' * hmat() #Rotation Matrix: [world_vector = Q * body_vector] 

    ## MAIN ROTOR FORCES
    Xmr = -u[3] * u[1] #Force along Body-X
    Ymr = u[3] * u[2] #Force along Body-Y
    Zmr = -u[3] #Force along Body-Z
    ## TAIL ROTOR FORCES
    Ytr = u[4]
    
    ## MAIN ROTOR MOMENTS
    Lmr = Ymr * hmr #Moment about Body-X
    Mmr = -Xmr * hmr #Moment about Body-Y
    ## MAIN ROTOR YAW TORQUE AGAINST AIR
    Qe = 0.0005 * rho * (omega_mr * Rmr)^2 * pi * Rmr^3  #(*Check Qe!)
    ## TAIL ROTOR MOMENTS
    Ntr = -Ytr * ltr #Moment about Body-Z
    Ltr = Ytr * htr #Moment about Body-X
    
    ## DRAG FORCES
    w_b = Q' * w #wind in Body Frame

    # Note: need to take care of the sign manually...
    Xfus = sign(w_b[1] - v[1]) * (0.5*rho*Sx_fus*Cd*(w_b[1] - v[1])*(w_b[1] - v[1])) #Drag force along Body-X
    Yfus = sign(w_b[2] - v[2]) * (0.5*rho*Sy_fus*Cd*(w_b[2] - v[2])*(w_b[2] - v[2])) #Drag force along Body-Y
    Zfus = sign(w_b[3] - v[3]) * (0.5*rho*Sz_fus*Cd*(w_b[3] - v[3])*(w_b[3] - v[3])) #Drag force along Body-Z

    ## GRAVITY IN BODY-FRAME
    g_matrix = Q' * [0; 0; -g]
    
    ## NEWTON-EULER'S EQUATIONS
    r_dot = Q * v
    q_dot = 0.5 * lmult(q) * hmat() * ω
    v_dot = SA[(v[2] * ω[3]) - (v[3] * ω[2]) - g_matrix[1] + ((Xmr + Xfus) / m);
             (v[3] * ω[1]) - (v[1] * ω[3]) + g_matrix[2] + ((Ymr + Yfus + Ytr) / m);
             (v[1] * ω[2]) - (v[2] * ω[1]) + g_matrix[3] + ((Zmr + Zfus) / m)]
    ω_dot = SA[((ω[2] * ω[3] * (Iyy - Izz)) + (Lmr + Ltr)) / Ixx;
             ((ω[1] * ω[3] * (Izz - Ixx)) + (Mmr)) / Iyy;
             ((ω[1] * ω[2] * (Ixx - Iyy)) + (-Qe + Ntr)) / Izz]


    vw_dot = vw #wind dynamics are not known...it does not matter

    return SA[r_dot; q_dot; v_dot; ω_dot; vw_dot; 0; 0; 0]  
    
end


"""
RobotDynamics.inertia(model::Helicopter) = model.J
RobotDynamics.inertia_inv(model::Helicopter) = model.Jinv
RobotDynamics.mass(model::Helicopter) = model.mass

function Base.zeros(model::Helicopter{R}) where R
    x = RobotDynamics.build_state(model, zero(RBState))
    u = @SVector fill()
    return x, u

end
"""



####################################################################
#### For Visualizer: ###############################################
####################################################################
####################################################################
"""
# Animation
function set_heli_model!(vis::Visualizer)
    #### Note: all dimensions in metres.
    width = 0.05
    # Length:
    Ltr = 0.91
    
    # Tail Rotor:
    htr = 0.08
    Rtr = 0.13
    
    # Main Rotor
    hmr = 0.235
    Rmr = 0.775
    
    # Fuselage:
    rad_fuse = 0.20
    Rshaft = 0.01

    fuselage = Sphere(Point3(0,0.0,0.0),rad_fuse)
    #fuselage2 = Sphere(Point3(-0.15,0.0,0.0),rad_fuse)
    MainRotor = Cylinder(Point3(0.0,0.0,hmr), Point3(0.0,0.0,hmr + 0.01), Rmr)
    Rotorshaft = Cylinder(Point3(0.0,0.0,0.0), Point3(0.0,0.0,hmr), Rshaft)
    Chassis = Rect3D(Vec(-Ltr,-width/2,-width/2), Vec(Ltr,width,width)) #body diagonal of a cuboid ??
    TailRotor = Cylinder(Point3(-Ltr,-0.04,htr), Point3(-Ltr,0.04,htr), Rtr)


    setobject!(vis["geom"]["fuselage"], fuselage, MeshPhongMaterial(color=colorant"gray"))
    #setobject!(vis["geom"]["fuselage2"], fuselage2, MeshPhongMaterial(color=colorant"gray"))
    setobject!(vis["geom"]["MainRotor"], MainRotor, MeshPhongMaterial(color=colorant"black"))
    setobject!(vis["geom"]["Rotorshaft"], Rotorshaft, MeshPhongMaterial(color=colorant"gray"))
    setobject!(vis["geom"]["Chassis"], Chassis, MeshPhongMaterial(color=colorant"blue"))
    setobject!(vis["geom"]["TailRotor"], TailRotor, MeshPhongMaterial(color=colorant"black"))
    
    return Nothing
end

function animate(vis, X, T, h)
    anim = MeshCat.Animation(Int(1/h))
    for t=1:size(thist,1) #10000
        atframe(anim, t) do
            #trans = Translation(10.0,0.0,0.0)
            trans = Translation(X[1,t],X[2,t],X[3,t])
            Q = qtoQ(X[4:7,t])
            rot = LinearMap(qtoQ([0.0,1.0,0.0,0.0]) * Q)
            settransform!(vis, trans ∘ rot)
        end
    end
    return anim
end

## Completed as a part of OCRL course project Spring-21: Sharvit Dabir, Vaibhav Shete, Advaith Sethuraman, Rishi Narayan Chakraborty.
## This project wouldn't have been possible without the help of Dr. Zachary Manchester and Brian Jackson.
"""