using VortexLattice
x = [0.5000000099993505, 0.7500000099998233, 0.2000000099996648, 0.4000000099991854, 0.29999999000018374, -9.998559696625576e-9, 0.3, 0.29999809265135313, 1.5999999840682646]
function tester(x)    
    function airframe_surface_initializer(wing_geo, horizontal_geo, vertical_geo, h_position, v_position)
        #generate panels for wing
        wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc,
        fc = fc_w, spacing_s = spacing_s, spacing_c = spacing_c, mirror = mirror)
    
        #generate panels for horizontal Stabilizer
        hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h,
        fc = fc_h, spacing_s = spacing_s_h, spacing_c = spacing_c_h, mirror = mirror_h)
        VortexLattice.translate!(hgrid, h_position)
        VortexLattice.translate!(htail, h_position)
    
        #generate panels for vertical Stabilizer
        vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v,
        fc = fc_v, spacing_s = spacing_s_v, spacing_c = spacing_c_v, mirror = mirror_v)
        VortexLattice.translate!(vgrid, v_position)
        VortexLattice.translate!(vtail, v_position)
    
        grids = [wgrid, hgrid, vgrid]
        surfaces = [wing, htail, vtail]
        
        return grids, surfaces
    end
    
    function lift_to_drag_solver(airframe_grid, airframe_surface, ref, fs, symmetric)
        surface_id = [1,2,3]
    
        #solve for coefficients using VLM
        system = steady_analysis(airframe_surface, ref, fs, symmetric = symmetric, surface_id = surface_id)
        CF, CM = body_forces(system, frame = Wind())
        CDiff = far_field_drag(system)
    
        CD, CY, CL = CF
        Cl, Cm, Cn = CM
    
        return CL/CD, CL
    end

    function stability_derivatives_solver(airframe_grid, airframe_surface, ref, fs, symmetric)
        surface_id = [1,2,3]
    
        #solve for derivatives using VLM
        system = steady_analysis(airframe_surface, ref, fs, symmetric = symmetric, surface_id = surface_id)
    
        dCF, dCM = stability_derivatives(system)
    
        CDa, CYa, CLa = dCF.alpha
        Cla, Cma, Cna = dCM.alpha
        CDb, CYb, CLb = dCF.beta
        Clb, Cmb, Cnb = dCM.beta
        CDp, CYp, CLp = dCF.p
        Clp, Cmp, Cnp = dCM.p
        CDq, CYq, CLq = dCF.q
        Clq, Cmq, Cnq = dCM.q
        CDr, CYr, CLr = dCF.r
        Clr, Cmr, Cnr = dCM.r
    
        return Cma, Cnb, Clb
    end
    
    x1, x2, x3, x4, x5, x6, x7, x8, x9 = x

    #Wing Parameters
    xle = [0.0, x1]
    yle = [0.0, x2]
    zle = [0.0, x3]
    chord = [0.3, 0.2]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc_w = fill((xc) -> 0, 2)
    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    wing_geo = (xle, yle, zle, chord, theta, phi, fc_w, ns, nc, spacing_s, spacing_c, mirror)

    #Horizontal Stabalizer Parameters
    xle_h = [0.0, x4]
    yle_h = [0.0, x5]
    zle_h = [0.0, x6]
    chord_h = [0.2, 0.2]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2)
    ns_h = 6
    nc_h = 3
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = false
    horizontal_geo = (xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h)

    #Vertical Stabilizer Parameters
    xle_v = [0.0, x7]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, x8]
    chord_v = [0.2, 0.1]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    fc_v = fill((xc) -> 0, 2)
    ns_v = 5
    nc_v = 3
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false
    vertical_geo = (xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v)

    #Reference Parameters
    Sref = 2*yle[2]*((chord[1]+chord[2])/2)
    cref = (chord[1]+chord[2])/2
    bref = 2*yle[2]
    rref = [0.50, 0.0, 0.0]
    Vinf = 22.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    #Freestream Parameters
    Vinf = 22.0
    alpha = 2*pi/180
    beta = 0.0
    Omega = [0.0, 0.0, 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    #Symmetry of Airframe Body
    symmetric = [true, true, false]

    #Horizontal and Vertical Stabilizer positioning
    h_position = [x9, 0.0, 0.0]
    v_position = [x9, 0.0, 0.0]

    airframe_grid, airframe_surface = airframe_surface_initializer(wing_geo, horizontal_geo, vertical_geo, h_position, v_position)

    lift_drag, lift = lift_to_drag_solver(airframe_grid, airframe_surface, ref, fs, symmetric)

    Cma, Cnb, Clb = stability_derivatives_solver(airframe_grid, airframe_surface, ref, fs, symmetric)


    
    return (Clb)

end

tester(x)