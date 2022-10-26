using VortexLattice
using Plots


"""
Pseudocode Discussion:

    Standard Steady State Analysis:
        First: Set up the geometry of the wing
        Second: Set discretization of the wing (e.g. how many panels and how these panels are spaced)
        Third: Set up reference parameters (e.g. Areas, span length, freestream velocity, etc.)
        Fourth: Set up freestream parameters (e.g. Angle of attack, side slip angle, etc.)
        Fifth: Use the method wing_to_surface_panels in order to generate the surface
        Sixth: Generate a vector with all the surface data (e.g. surface = [surface])
        Seventh: If the flow conditions are symmetric about the X-Z axis, We can declare that symmetric = true
        Eigth: Use the method steady_analysis in order to start finding coefficients
        Ninth: Perform a far field ananlysis using far_field_drag in order to get a better value for the induced drag coefficient
    For Asymmetric Flow Conditions:
        *** Everything is basically the same as the standard steady state analysis except that we set mirror = true in the
        wing_to_surface_panels method.
        *** Also, make sure to set symmetric to false.
    Finding Stability Derivatives:
        *** This is not too hard to find after calculating CF and CM; we use the stability_derivatives method. 
    Steady State for a Dihedral Wing:
        *** It is the same set up as the standard steady state except that when setting up the geometry, give the Z-axis length a value other than 
        0.
    Steady state for Wing and Tail:
        *** It is the same set up as the standard steady state, but we add the geometry and discretization for tail, as well. 
        *** Normal vectors for the new geometries need to be set up, too? 
    Acceleration Simulation:
        --- The idea makes sense, but the application and method are interesting ---
"""

#Wing Geometry
#Set Up Geometry
xle = [0.0, 1.0]
yle = [0.0,7.5]
zle = [0.0,2.0]
chord = [2.2, 1.0]
theta = [0.0,0.0]
phi = [0.0,0.0]
fc = fill((xc) -> 0, 2)
#Discretization Values
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false
Wing_Geo = (xle, yle, zle, chord, theta, phi, ns, nc, spacing_s, spacing_c, mirror, fc)

#Horizontal Stabalizer Geometry
xle_h = [0.0,0.14]
yle_h = [0.0,1.25]
zle_h = [0.0,0.0]
chord_h = [0.7,0.42]
theta_h = [0.0,0.0]
phi_h = [0.0,0.0]
fc_h = fill((xc) -> 0, 2)
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = false
H_Tail_Geo = (xle_h, zle_h, yle_h, chord_h, theta_h, phi_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, fc_h)

#Vertical Stabalizer Geometry
xle_v = [0.0,0.14]
yle_v = [0.0,0.0]
zle_v = [0.0,1.0]
chord_v = [0.7, 0.42]
theta_v = [0.0,0.0]
phi_v = [0.0,0.0]
fc_v = fill((xc) -> 0, 2)
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false
V_Tail_Geo = (xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, fc_v)



#Reference Parameters

Sref = 9.0
cref = 0.9
bref = 10.0
rref = [0.50, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref,cref,bref,rref,Vinf)


#Freestream Parameters
Vinf = 1.0
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0,0.0,0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = [true,true,false]

"""
    VLM_Solver_Basic(xle, yle, zle, chord, theta, phi, ns, nc, fc, spacing_s, spacing_c, symmetric, ref)

Finds the aerodynamic coefficients with an input of just a wing.

The parameters of this function are just single inputs instead of the tuples that are found in the other functions.

# Example of input parameters:

xle_h = [0.0,0.14] 'This is the distance that the leading edge is pushed back'

yle_h = [0.0,1.25] 'This is the span of a single wing'

zle_h = [0.0,0.0] 'This is how high in the dihedral the direction the wing goes'

chord_h = [0.7,0.42] 'This is the length of the root chord to the tip chord'

theta_h = [0.0,0.0] 'This is the rotation of the wing in an angle of attack manner'

phi_h = [0.0,0.0] 'This is the twist along the wing'

fc_h = fill((xc) -> 0, 2) 'This is the camber along the wing'

ns_h = 6 'This is how many sections/panels are created in the spanwidth direction'

nc_h = 3 'This is how many sections/panels are created in the chordwidth direction'

spacing_s_h = Uniform() 'This describes the spacing of the created panels in the spanwidth direction'

spacing_c_h = Uniform() 'This describes the spacing of the created panels in the chordwidth direction'

mirror_h = false 

symmetric = true 'This describes whether the flow is identical on either side of the aircraft'

ref = Reference(Sref,cref,bref,rref,Vinf)
"""
function vlm_solver_basic(xle, yle, zle, chord, theta, phi, ns, nc, fc, spacing_s, spacing_c, symmetric, ref)
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc,
    fc = fc, spacing_s = spacing_s, spacing_c = spacing_c)
    surfaces = [surface]

    #Perform steady state analysis
    system = steady_analysis(surfaces, ref, fs, symmetric = symmetric)

    #Retrieve near field forces
    CF, CM = body_forces(system, frame=Wind())
    
    #Perform Far Field Analysis
    CDiff = far_field_drag(system)

    CD,CY,CL = CF
    Cl,Cm,Cn = CM

    println("Coefficient of Induced Drag = " ,CD)
    println("Coefficient of Lift = " ,CL)
end


"""
    VLM_Solver_Efficiency(xle, zle, chord, theta, phi, ns, nc, fc, spacing_s, spacing_c, symmetric, fs)

Similar to the VLM_Solver_Basic but computes the inviscid efficiency as well as the coefficients.

The parameters of this function are just single inputs instead of the tuples that are found in the other functions.

# Example of input parameters:
    
xle_h = [0.0,0.14] 'This is the distance that the leading edge is pushed back'

yle_h = [0.0,1.25] 'This is the span of a single wing'

zle_h = [0.0,0.0] 'This is how high in the dihedral the direction the wing goes'

chord_h = [0.7,0.42] 'This is the length of the root chord to the tip chord'

theta_h = [0.0,0.0] 'This is the rotation of the wing in an angle of attack manner'

phi_h = [0.0,0.0] 'This is the twist along the wing'

fc_h = fill((xc) -> 0, 2) 'This is the camber along the wing'

ns_h = 6 'This is how many sections/panels are created in the spanwidth direction'

nc_h = 3 'This is how many sections/panels are created in the chordwidth direction'

spacing_s_h = Uniform() 'This describes the spacing of the created panels in the spanwidth direction'

spacing_c_h = Uniform() 'This describes the spacing of the created panels in the chordwidth direction'

mirror_h = false 

symmetric = true 'This describes whether the flow is identical on either side of the aircraft'

ref = Reference(Sref,cref,bref,rref,Vinf)
"""
function vlm_solver_efficiency(xle, zle, chord, theta, phi, ns, nc, fc, spacing_s, spacing_c, symmetric, fs)
    Aspect_Ratios = []
    efficiencies = []
    for i in 10:0.5:1000
        Sref = 2*i*chord[1]
        cref = chord[1]
        bref = 2*i
        rref = [0.50, 0.0, 0.0]
        Vinf = 1.0
        ref = Reference(Sref,cref,bref,rref,Vinf)
        yle = [0.0 , i]
        grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc,
        fc = fc, spacing_s = spacing_s, spacing_c = spacing_c)
        surfaces = [surface]

        #Perform steady state analysis
        system = steady_analysis(surfaces, ref, fs, symmetric = symmetric)

        #Retrieve near field forces
        CF, CM = body_forces(system, frame=Wind())
        
        #Perform Far Field Analysis
        CDiff = far_field_drag(system)

        CD,CY,CL = CF
        Cl,Cm,Cn = CM

        #Span length and Reference Area
        b = bref
        S = Sref


        #Aspect Ratio
        AR = (b^2)/S


        #efficiency equation
        e = ((CL)^(2))/(pi*AR*CD)
        push!(efficiencies , e)
        push!(Aspect_Ratios, AR)
    end
    plot(Aspect_Ratios, efficiencies, xlabel = "Aspect Ratio", ylabel = "Inviscid Efficiency")
end

"""
    VLM_Solver_Wing_Tail_Derivatives(Wing_Geo, H_Tail_Geo, V_Tail_Geo, ref, fs, symmetric)

This function implements a tail as well as a wing into the geometry being analyzed. Additionally, this function calculates the stability derivatives
of the airframe.

The parameters of this function are inputed as tuples instead of the single inputs put into the wing only functions.

# Example of inputed tuple

xle_h = [0.0,0.14] 'This is the distance that the leading edge is pushed back'

yle_h = [0.0,1.25] 'This is the span of a single wing'

zle_h = [0.0,0.0] 'This is how high in the dihedral the direction the wing goes'

chord_h = [0.7,0.42] 'This is the length of the root chord to the tip chord'

theta_h = [0.0,0.0] 'This is the rotation of the wing in an angle of attack manner'

phi_h = [0.0,0.0] 'This is the twist along the wing'

fc_h = fill((xc) -> 0, 2) 'This is the camber along the wing'

ns_h = 6 'This is how many sections/panels are created in the spanwidth direction'

nc_h = 3 'This is how many sections/panels are created in the chordwidth direction'

spacing_s_h = Uniform() 'This describes the spacing of the created panels in the spanwidth direction'

spacing_c_h = Uniform() 'This describes the spacing of the created panels in the chordwidth direction'

mirror_h = false 

H_Tail_Geo = (xle_h, zle_h, yle_h, chord_h, theta_h, phi_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, fc_h)
"""
function vlm_solver_wing_tail_derivatives(Wing_Geo, H_Tail_Geo, V_Tail_Geo, ref, fs, symmetric)
    #generate panels for wing
    wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc,
    fc = fc, spacing_s = spacing_s, spacing_c = spacing_c, mirror = mirror)
        
    #generate panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h,
    fc = fc_h, spacing_s = spacing_s_h, spacing_c = spacing_c_h, mirror = mirror_h)
    VortexLattice.translate!(hgrid, [4.0,0.0,0.0])
    VortexLattice.translate!(htail, [4.0,0.0,0.0])

    #generate panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v,
    fc = fc_v, spacing_s = spacing_s_v, spacing_c = spacing_c_v, mirror = mirror_v)
    VortexLattice.translate!(vgrid, [4.0,0.0,0.0])
    VortexLattice.translate!(vtail, [4.0,0.0,0.0])

    grids = [wgrid, hgrid, vgrid]
    surfaces = [wing, htail, vtail]

    surface_id = [1,2,3]

    system = steady_analysis(surfaces, ref, fs, symmetric = symmetric, surface_id = surface_id)
    CF , CM = body_forces(system, frame = Wind())
    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    #Finding the Stability Derivatives

    dCF, dCM = stability_derivatives(system)
end


"""
    VLM_Solver_Wing_Tail(Wing_Geo, H_Tail_Geo, V_Tail_Geo, ref, symmetric)

This function implements a tail as well as a wing into the geometry being analyzed. This function can be used to find the aerodynamic coefficients, but it
doesn't solve for the stability derivatives.

The parameters of this function are inputed as tuples instead of the single inputs put into the wing only functions.

# Example of inputed tuple

xle_h = [0.0,0.14] 'This is the distance that the leading edge is pushed back'

yle_h = [0.0,1.25] 'This is the span of a single wing'

zle_h = [0.0,0.0] 'This is how high in the dihedral the direction the wing goes'

chord_h = [0.7,0.42] 'This is the length of the root chord to the tip chord'

theta_h = [0.0,0.0] 'This is the rotation of the wing in an angle of attack manner'

phi_h = [0.0,0.0] 'This is the twist along the wing'

fc_h = fill((xc) -> 0, 2) 'This is the camber along the wing'

ns_h = 6 'This is how many sections/panels are created in the spanwidth direction'

nc_h = 3 'This is how many sections/panels are created in the chordwidth direction'

spacing_s_h = Uniform() 'This describes the spacing of the created panels in the spanwidth direction'

spacing_c_h = Uniform() 'This describes the spacing of the created panels in the chordwidth direction'

mirror_h = false 
    
H_Tail_Geo = (xle_h, zle_h, yle_h, chord_h, theta_h, phi_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, fc_h)
"""
function vlm_solver_wing_tail(Wing_Geo, H_Tail_Geo, V_Tail_Geo, ref, symmetric)
    alpha_list = []
    Lift_C = []
    for i in -15.0:15.0
        alpha = (i*pi/180)
        Vinf = 1.0
        beta = 0.0
        Omega = [0.0,0.0,0.0]
        fs = Freestream(Vinf, alpha, beta,Omega)
        #generate panels for wing
        wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc,
        fc = fc, spacing_s = spacing_s, spacing_c = spacing_c, mirror = mirror)
            
        #generate panels for horizontal tail
        hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h,
        fc = fc_h, spacing_s = spacing_s_h, spacing_c = spacing_c_h, mirror = mirror_h)
        VortexLattice.translate!(hgrid, [4.0,0.0,0.0])
        VortexLattice.translate!(htail, [4.0,0.0,0.0])

        #generate panels for vertical tail
        vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v,
        fc = fc_v, spacing_s = spacing_s_v, spacing_c = spacing_c_v, mirror = mirror_v)
        VortexLattice.translate!(vgrid, [4.0,0.0,0.0])
        VortexLattice.translate!(vtail, [4.0,0.0,0.0])

        grids = [wgrid, hgrid, vgrid]
        surfaces = [wing, htail, vtail]

        surface_id = [1,2,3]

        system = steady_analysis(surfaces, ref, fs, symmetric = symmetric, surface_id = surface_id)
        CF , CM = body_forces(system, frame = Wind())
        CDiff = far_field_drag(system)

        CD, CY, CL = CF
        Cl, Cm, Cn = CM
        push!(alpha_list, alpha)
        push!(Lift_C, CL)
    end
    plot(alpha_list, Lift_C, xlabel = "Angle of Attack In Degrees (Body of the Plane)", ylabel = "Coefficient of Lift")
end

