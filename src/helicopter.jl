########################################################################
#### CONTENTS:##########################################################
## hat(v)
## L(q)
## qtoQ(q)
## helicopter_dynamics(x, u)
## set_heli_model!(vis::Visualizer) **commented out
## animate(vis, X, T, h) **commented out
########################################################################


########################################################################
#### Quaternion Functions: #############################################
########################################################################
########################################################################
function hat(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function L(q)
    ## Computes Left Multiplier of Quaternion:
    s = q[1]
    v = q[2:4]
    L = [s -v'; v (s * I + hat(v))]
    return L
end

function qtoQ(q)
    ## Converts Quaternion to Rotation Matrix
    T = Diagonal([1; -ones(3)])
    H = [zeros(1,3); I]
    return H'*T*L(q)*T*L(q)*H
end

####################################################################
#### X-CELL 60 SE Helicopter Dynamics: #############################
####################################################################
####################################################################
function helicopter_dynamics(x, u)
    # x: [x y z Quat u v w p q r wx wy wz wxdot wydot wzdot]'
    # u: [a b Tmr Ttr]'

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
    ## GRAVITY
    g = 9.81 #m/s^2

    ## MASS
    M = 8.2 #kg

    ## INERTIA
    Ixx = 0.18 #kg m^2
    Iyy = 0.34 #kg m^2
    Izz = 0.28 #kg m^2
    J = Diagonal([Ixx, Iyy, Izz])

    ## MAIN ROTOR
    hmr = 0.235 #m
    Ωmr = 167 #rad/s
    Rmr = 0.775 #m

    ## DRAG PARAMETERS
    # https://www.engineeringtoolbox.com/drag-coefficient-d_627.html
    Cd = 1.17 #square flat plate 90degrees to the wind.

    Sx_fus = 0.1 #m^2
    Sy_fus = 0.22 #m^2
    Sz_fus = 0.15 #m^2

    ## TAIL ROTOR
    ntr = 4.66 #Gear ratio of tr to mr
    Ωtr = ntr * Ωmr #rad/s
    Rtr = 0.13 #m
    ltr = 0.91 #m
    htr = 0.08 #m

    ## AIR PARAMETERS
    rho = 1.204 #kg/m^3 #ρ
    

    ####################################################################
    #### X-CELL 60 SE Helicopter Dynamics: #############################
    ####################################################################
    ####################################################################
    #=   
    #### NOTE ON STATES::                                                   
    x[14:16] --> wind velocities (estimated or ground-truth) 
    x[14:16] --> wind "accelerations" 

    #### NOTE ON DYNAMICS:
    ASSUME: No blade-flapping dynamics;
            No horizontal+Vertical fins;
            No stabilizer;
            Centre of Pressure == C.O.G

    #### NOTE ON CONTROL INPUTS:
    u[1] = Main Rotor tilt angle with vertical in Helicopter's X-Z plane (*check!) 
    u[2] = Main Rotor tilt angle with vertical in Helicopter's Y-Z plane
    u[3] = Main Rotor Thrust
    u[4] = Tail Rotor Thrust
    =#
    ####################################################################
    r = x[1:3] #x y z
    q = x[4:7] # quaternion
    v = x[8:10] #u v w
    ω = x[11:13] #p q r
    w = x[14:16] #wind-velocities
    vw = x[17:19] #wind-velocity_dot

    Q = qtoQ(q) #Rotation Matrix: [world_vector = Q * body_vector]
    
    ## MAIN ROTOR FORCES
    Xmr = -u[3] * u[1] #Force along Body-X
    Ymr = u[3] * u[2] #Force along Body-Y
    Zmr = -u[3] #Force along Body-Z
    
    ## MAIN ROTOR MOMENTS
    Lmr = Ymr * hmr #Moment about Body-X
    Mmr = -Xmr * hmr #Moment about Body-Y

    ## MAIN ROTOR YAW TORQUE AGAINST AIR
    Qe = 0.0005 * rho * (Ωmr * Rmr)^2 * pi * Rmr^3  #(*Check Qe!)
    
    ## TAIL ROTOR FORCES
    Ytr = u[4]
    
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
    #ṙ, q̇, v̇, ω̇ 
    r_dot = Q * v
    q_dot = 0.5 * L(q) * H * ω
    v_dot = [(v[2] * ω[3]) - (v[3] * ω[2]) - g_matrix[1] + ((Xmr + Xfus) / M);
             (v[3] * ω[1]) - (v[1] * ω[3]) + g_matrix[2] + ((Ymr + Yfus + Ytr) / M);
             (v[1] * ω[2]) - (v[2] * ω[1]) + g_matrix[3] + ((Zmr + Zfus) / M)]
    
    ω_dot = [((ω[2] * ω[3] * (Iyy - Izz)) + (Lmr + Ltr)) / Ixx;
             ((ω[1] * ω[3] * (Izz - Ixx)) + (Mmr)) / Iyy;
             ((ω[1] * ω[2] * (Ixx - Iyy)) + (-Qe + Ntr)) / Izz]
    
    #return [r_dot; q_dot; v_dot; ω_dot; vw; 0; 0; 0]    
end


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