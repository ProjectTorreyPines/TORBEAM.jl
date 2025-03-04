module TORBEAM
using IMAS
using Interpolations

    @Base.kwdef struct TorbeamParams
        npow = 1
        ncd = 1
        ncdroutine = 2
        nprofv = 50
        noout = 0
        nrela = 1
        nmaxh = 3
        nabsroutine = 1
        nastra = 0
        nprofcalc = 1
        ncdharm = 1
        nrel = 0

        #Float parameters
        xrtol = 1e-07
        xatol = 1e-07
        xstep = 2.0
        rhostop = 0.96
        xzsrch = 0
        
    end


    function torbeam!(dd::IMAS.dd, torbeam_params::TorbeamParams)
        # TORBEAM subroutine in IMAS
        #----------------------------
        # TORBEAM subroutine in IMAS
        #----------------------------
        
        # Define parameters
                # from libtorbeam/src/libsrc/dimensions.f90 -> Could be moved to a torbeam.yaml (?)
        # Define sizes of various arrays. Need to allocate them here to then pass to TORBEAM

        nbeam = length(dd.ec_launchers.beam)
        if nbeam < 1
            println("No beams - returning")
            return
        end
        maxint = 50
        maxflt = 50
        maxrhr = 20
        ndat = 100000
        npnt = 5000
        mmax = 450
        maxvol = 100
        ntraj = 10000
        nmax = mmax
        maxdim = 1+mmax+nmax+4*mmax*nmax
        maxlen = 2*mmax+2*nmax

        # Allocate input arrays
        rhoresult = zeros(Float64, maxrhr)
        intinbeam = zeros(Int32, maxint)
        floatinbeam = zeros(Float64, maxflt)
        eqdata = zeros(Float64, maxdim)
        prdata = zeros(Float64, maxlen)
        volprof = zeros(Float64, 2*maxvol)

        # Define output scalars
        iend = Ref{Int32}(0)
        kend = Ref{Int32}(0)
        icnt = Ref{Int32}(0)
        ibgout = Ref{Int32}(0)
        
        

        # Allocate output arrays
        t1data = zeros(Float64, 6*ndat)
        t1tdata = zeros(Float64, 6*ndat)
        t2data = zeros(Float64, 5*ndat)
        t2ndata = zeros(Float64, 3*npnt)
        trajout = zeros(Float64, (nbeam, 15, ntraj))

        # Allocate arrays
        eqt2d = dd.equilibrium
        eqt2d = IMAS.findfirst(:rectangular, eqt2d.time_slice[].profiles_2d)
        psi2d = eqt2d.psi
        Rarr = eqt2d.grid.dim1
        Zarr = eqt2d.grid.dim2
        ni = length(Rarr)
        nj = length(Zarr)
        br = eqt2d.b_field_r
        bt = eqt2d.b_field_tor
        bz = eqt2d.b_field_z

        # Interpolator for psi -> rho-tor need that later
        eq1d = dd.equilibrium.time_slice[].profiles_1d
        rho_tor_interpolator = Interpolations.linear_interpolation(eq1d.psi, eq1d.rho_tor_norm)

        psiedge = dd.equilibrium.time_slice[].global_quantities.psi_boundary
        psiax = dd.equilibrium.time_slice[].global_quantities.psi_axis
        cp1d = dd.core_profiles.profiles_1d[]
        # Initialize vector 'eqdata' as TORBEAM input (topfile):
        # Psi, 1D, for normalization
        npsi = length(cp1d.grid.psi)
        psi = cp1d.grid.psi
        # To ensure consistency between the 1D and 2D psi profiles: take both from the equilibrium IDS

        # FILL TORBEAM INTERNAL EQUILIBRIUM DATA
        eqdata[1] = psiedge
        for i = 1:ni
            eqdata[i+1] = Rarr[i]
        end
        for j = 1:nj
            eqdata[j+ni+1] = Zarr[j]
        end
        k = 0
        for j = 1:nj
            for i = 1:ni
                k = k+1
                eqdata[k+ni+nj+1] = br[i,j]
            end
        end
        k = 0
        for j = 1:nj
            for i = 1:ni
                k = k+1
                eqdata[k+ni+nj+ni*nj+1] = bt[i,j]
            end
        end
        k = 0
        for j = 1:nj
            for i = 1:ni
                k = k+1
                eqdata[k+ni+nj+2*ni*nj+1] = bz[i,j]
            end
        end
        k = 0
        for j = 1:nj
            for i = 1:ni
                k = k+1
                eqdata[k+ni+nj+3*ni*nj+1] = psi2d[i,j]
            end
        end

        iplasma = dd.equilibrium.time_slice[].global_quantities.ip

        # FILL TORBEAM INTERNAL PROFILE DATA
        #... Initialize vector 'prdata' as TORBEAM input (ne.dat & Te.dat):
        # Psi and profiles
        te = dd.core_profiles.profiles_1d[].electrons.temperature
        ne = dd.core_profiles.profiles_1d[].electrons.density

        for i = 1:npsi
            prdata[i] = sqrt((psi[i]-psi[1])/(psi[npsi]-psi[1]))
            prdata[i+npsi] = ne[i]*1.0e-19
        end
        for i = 1:npsi
            prdata[i+2*npsi] = sqrt((psi[i]-psi[1])/(psi[npsi]-psi[1]))
            prdata[i+2*npsi+npsi] = te[i]*1.0e-3
        end
        # Initialize antenna data as TORBEAM input:

        # DETERMINE WHETHER PSI FLUX IS MAXIMUM (1) OR MINIMUM (-1) AT THE MAGNETIC AXIS
        if (psiedge > psiax)
            sgnm = 1.
        else
            sgnm = -1.
        end
        npointsout = zeros(Int64, nbeam)
        extrascal = zeros(Float64, (nbeam,5))
        profout = zeros(Float64, (nbeam, 3, npnt))
        # LOOP OVER BEAMS OF THE EC_LAUNCHERS IDS
        for ibeam = 1:nbeam

            # ONLY DEAL WITH ACTIVE BEAMS
            if (@ddtime(dd.ec_launchers.beam[ibeam].power_launched.data) > 0)

                # IT LOOKS LIKE TORBEAM NEEDS PHI = 0, OTHERWISE IT DOES NOT TREAT THE BEAM PROPERLY
                # BUT WE WILL RESTORE THE ACTUAL PHI ANGLE AFTER THE RAY-TRACING, SO WE DON'T
                # PUT ec_launchers%BEAM(IBEAM)%LAUNCHING_POSITION%PHI TO 0 ANYMORE
                # (WE ARTIFICIALLY PUT PHI=0 IN floatinbeam(3) AND FLOTINBEAM(4) INSTEAD

                #intinbeam
                intinbeam[1] = 2  # tbr
                intinbeam[2] = 2  # tbr
                intinbeam[3] = dd.ec_launchers.beam[ibeam].mode  # (nmod)
                intinbeam[4] = torbeam_params.npow
                intinbeam[5] = torbeam_params.ncd
                intinbeam[6] = 2  # tbr
                intinbeam[7] = torbeam_params.ncdroutine
                intinbeam[8] = torbeam_params.nprofv
                intinbeam[9] = torbeam_params.noout
                intinbeam[10] = torbeam_params.nrela
                intinbeam[11] = torbeam_params.nmaxh
                intinbeam[12] = torbeam_params.nabsroutine
                intinbeam[13] = torbeam_params.nastra
                intinbeam[14] = torbeam_params.nprofcalc
                intinbeam[15] = torbeam_params.ncdharm
                intinbeam[16] = 0
                intinbeam[17] = 0
                intinbeam[maxint] = torbeam_params.nrel

                #floatinbeam(17:18): obsolete --> not filled)
                #floatinbeam(6:13):  analytic --> not filled)
                #floatinbeam(26:32): analytic --> not filled)
                floatinbeam[1] = @ddtime(dd.ec_launchers.beam[ibeam].frequency.data)  # (xf)
                # floatinbeam[2] = rad2deg(-@ddtime(dd.ec_launchers.beam[ibeam].steering_angle_tor))
                # floatinbeam[3] = rad2deg(@ddtime(dd.ec_launchers.beam[ibeam].steering_angle_pol))
                phi_tor  = @ddtime(dd.ec_launchers.beam[ibeam].steering_angle_tor)
                theta_pol =  @ddtime(dd.ec_launchers.beam[ibeam].steering_angle_pol)
                floatinbeam[2]  = -np.rad2deg(np.arcsin(np.cos(theta_pol)*np.sin(phi_tor)))
                floatinbeam[3]  = np.rad2deg(np.arctan2(np.tan(theta_pol), np.cos(phi_tor)))
                floatinbeam[4] = 1.e2*dd.ec_launchers.beam[ibeam].launching_position.r[1] * cos(0)  # (xxb)
                floatinbeam[5] = 1.e2*dd.ec_launchers.beam[ibeam].launching_position.r[1] * sin(0)  # (xyb)
                floatinbeam[6] = 1.e2*dd.ec_launchers.beam[ibeam].launching_position.z[1]  # (xzb)

                floatinbeam[15] = torbeam_params.xrtol  # keep
                floatinbeam[16] = torbeam_params.xatol  # keep
                floatinbeam[17] = torbeam_params.xstep  # keep
                floatinbeam[20] = -1.e2/(dd.ec_launchers.beam[ibeam].phase.curvature[1,1])  # (xryyb)
                floatinbeam[21] = -1.e2/(dd.ec_launchers.beam[ibeam].phase.curvature[2,1])  # (xrzzb)
                if (cos(@ddtime(dd.ec_launchers.beam[ibeam].spot.angle))^2 > 0.5)
                    floatinbeam[22] = dd.ec_launchers.beam[ibeam].spot.size[1,1]*1.e2  # (xwyyb)
                    floatinbeam[23] = dd.ec_launchers.beam[ibeam].spot.size[2,1]*1.e2  # (xwzzb)
                else
                    floatinbeam[22] = dd.ec_launchers.beam[ibeam].spot.size[2,1]*1.e2  # (xwzzb)
                    floatinbeam[23] = dd.ec_launchers.beam[ibeam].spot.size[1,1]*1.e2  # (xwyyb)
                end
                floatinbeam[24] = @ddtime(dd.ec_launchers.beam[ibeam].power_launched.data)*1.e-6  # (xpw0)
                floatinbeam[25] = dd.equilibrium.time_slice[].boundary.geometric_axis.r*1e2  # (xrmaj)
                floatinbeam[26] = dd.equilibrium.time_slice[].boundary.minor_radius*1e2  # (xrmin)
                floatinbeam[27] = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)
                floatinbeam[34] = sgnm  # (deduced from psi_ed-psi_ax)
                floatinbeam[35] = dd.core_profiles.profiles_1d[].zeff[1]  # (xzeff)
                floatinbeam[36] = torbeam_params.rhostop  # keep
                floatinbeam[37] = torbeam_params.xzsrch  # keep

                println("------------------------------------------------------------")
                println("Input power beam: ", ibeam, " ", @ddtime(dd.ec_launchers.beam[ibeam].power_launched.data)*1.e-6, " MW")

                # CALL TORBEAM
                ccall(
                    (:beam_, get(ENV, "TORBEAM_DIR", "") * "/../lib/libtorbeamC1.so"),    # Name in the shared library (append `_`)
                    Cvoid,                             # Return type
                    (Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, # Inputs
                     Ptr{Float64}, Ref{Cint}, Ptr{Float64}, Ptr{Float64}, Ref{Cint},
                     Ptr{Float64}, Ptr{Float64}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{Float64}), # Argument types
                    intinbeam, floatinbeam, ni, nj, eqdata, npsi, npsi, prdata,
                    rhoresult, iend, t1data, t1tdata, kend,
                    t2data, t2ndata, icnt, ibgout, torbeam_params.nprofv, volprof
                )
                println(rhoresult[19])
                # --------------------------------------------------------------------------------------------------
                # Result structures:
                # --------------------------------------------------------------------------------------------------
                # Beam propagation:
                # - t1data  = (6 variables)
                #   * R - major-radius coordinate of the central ray                  (0:iend-1)
                #   * Z - vertical coordinate of the central ray                      (iend:2*iend-1)
                #   * R - major radius of the upper peripheral ray, i.e.interaction   (2*iend:3*iend-1)
                #         of beam width with the poloidal plane above the central ray
                #   * Z - vertical coordinate of the upper peripheral ray             (3*iend:4*iend-1)
                #   * R - major radius of the lower peripheral ray                    (4*iend:5*iend-1)
                #   * Z - vertical coordinate of the lower peripheral ray             (5*iend:6*iend-1).
                # --------------------------------------------------------------------------------------------------
                # - t1tdata = (4 variables)
                #   * X-coordinate of the central ray
                #   * Y-coordinate of the central ray (i.e. projection of the central ray onto a horizontal plane)
                #   * X-coordinate of left and right peripheral rays
                #   * Y-coordinate of left and right peripheral rays
                #   (intersection of the beam width and the horizontal plane running throuh the central ray)
                # --------------------------------------------------------------------------------------------------
                # Absorption and current drive
                # --------------------------------------------------------------------------------------------------
                # - t2data structure (nprofv = number of radial points in profiles)
                # The first 3*nprofv entries of t2data are already taken by the radial profile
                # of the area of the flux surfaces and their volume
                # --> Other variables (scalars) start from rhoresult(4) = real(3*nprofv)
                # * GROUP-VELOCITY COMPONENTS:
                #   (dimensionless components of a unit vector tangent to the propagation direction of the central ray)
                #   t2data(3*nprofv)   = vx/denth (central ray)
                #   t2data(3*nprofv+1) = vy/denth (idem)
                #   t2data(3*nprofv+2) = vz/denth (idem)
                # * WIDTHS AND CURVATURES:
                #   (wyb/wzb = distance between the central ray and the peripheral "rays" in the horiz/vert. direction)
                #   (1/syb, 1/szb = radii of curvature in the horizontal/vertical direction, which defines how far
                #    ahead - or behind, depending on the sign - the corresponding geometrical-optics ray would cross)
                #   t2data(3*nprofv+3) = wyb (from central ray to others)
                #   t2data(3*nprofv+4) = wzb (idem)
                #   t2data(3*nprofv+5) = 1.e0_rkind/syb (central ray: distance between actual and geometric ray)
                #   t2data(3*nprofv+6) = 1.e0_rkind/szb (idem)
                # * PRINCIPAL WIDTHS (for a check of the area):
                #   t2data(3*nprofv+7) = wmaj
                #   t2data(3*nprofv+8) = wmin
                # ---------------------------------------------------------------------------------
                # - t2ndata = (3 variables)
                #   * radial coordinate (rho_p or rho_t as above, 0:npnt-1)
                #   * power density in MW/m3 (npnt:2*npnt-1)
                #   * driven current density in MA/m2 (2*npnt:3*npnt-1)
                # --------------------------------------------------------------------------------------------------

                npointsout[ibeam] = iend[]
                extrascal[ibeam, 1] = dd.ec_launchers.beam[ibeam].power_launched.data[1]*1.e-6
                extrascal[ibeam, 2] = 1.e6*rhoresult[13]
                extrascal[ibeam, 3] = 1.e6*rhoresult[12]


                # 1D PROFILES OF RHO, DP/DV, J (CHECKED OK, DIM = NPNT = 5000)
                for kp = 0:2
                    for lfd = 1:npnt
                        profout[ibeam, kp + 1, lfd] = t2ndata[kp*npnt+lfd]
                    end
                end

                # TO PREVENT GOING BEYOND THE PRE-DEFINED TRAJOUT ARRAY
                if (iend[] > ntraj)
                    iend[] = ntraj
                end

                # TRAJECTORY OF CENTRAL RAY AND 4 PERIPHERAL "RAYS"
                for lfd = 1:iend[]
                    # 1st ray
                    trajout[ibeam, 1, lfd] = t1data[lfd]                    # r
                    trajout[ibeam, 2, lfd] = t1data[iend[]+lfd]               # z
                    #trajout[ibeam, 3, lfd] = acos(t1tdata[lfd]/t1data[lfd]) # phi
                    trajout[ibeam, 3, lfd] = atan(t1tdata[iend[]+lfd],t1data[lfd]) # phi
                    # 2nd ray
                    trajout[ibeam, 4, lfd] = t1data[2*iend[]+lfd]
                    trajout[ibeam, 5, lfd] = t1data[3*iend[]+lfd]
                    #trajout[ibeam, 6, lfd] = acos(t1tdata[lfd]/t1data[2*iend[]+lfd])
                    trajout[ibeam, 6, lfd] = atan(t1tdata[iend[]+lfd],t1data[lfd])
                    # 3rd ray
                    trajout[ibeam, 7, lfd] = t1data[4*iend[]+lfd]
                    trajout[ibeam, 8, lfd] = t1data[5*iend[]+lfd]
                    #trajout[ibeam, 9, lfd] = acos(t1tdata[lfd]/t1data[4*iend[]+lfd])
                    trajout[ibeam, 9, lfd] = atan(t1tdata[iend[]+lfd],t1data[iend[]+lfd])
                    # 4th ray
                    trajout[ibeam, 10, lfd] = sqrt(t1tdata[2*iend[]+lfd]^2+t1tdata[3*iend[]+lfd]^2)
                    trajout[ibeam, 11, lfd] = t1data[iend[]+lfd]
                    #trajout[ibeam, 12, lfd] = acos(t1tdata[2*iend[]+lfd]/trajout[ibeam, 10, lfd])
                    trajout[ibeam, 12, lfd] = atan(t1tdata[3*iend[]+lfd],t1tdata[2*iend[]+lfd])
                    # 5th ray
                    trajout[ibeam, 13, lfd] = sqrt(t1tdata[4*iend[]+lfd]^2+t1tdata[5*iend[]+lfd]^2)
                    trajout[ibeam, 14, lfd] = t1data[iend[]+lfd]
                    #trajout[ibeam, 15, lfd] = acos(t1tdata[4*iend[]+lfd]/trajout[ibeam, 13, lfd+1])
                    trajout[ibeam, 15, lfd] = atan(t1tdata[5*iend[]+lfd],t1tdata[4*iend[]+lfd])
                end

            end # TEST BEAM_POWER > 0

        end # LOOP OVER BEAMS OF EC_LAUNCHERS IDS

        # ----------------------------
        # SAVE RESULTS INTO WAVES IDS
        # ----------------------------

        # GLOBAL QUANTITIES
        dd.waves.ids_properties.homogeneous_time = 1
        dd.waves.time = zeros(Float64,1)
        dd.waves.time[1] = @ddtime(dd.equilibrium.time)
        dd.waves.code.name = "TORBEAM"

        source_index = length(dd.core_sources.source) + 1
        resize!(dd.core_sources.source, source_index)
        dd.core_sources.source[source_index].identifier.index = 3
        resize!(dd.core_sources.source[source_index].profiles_1d, 1)
        dd.core_sources.source[source_index].profiles_1d[1].time = @ddtime(dd.equilibrium.time)

        # LOOP OVER BEAMS (LAUNCHERS)
        #if(nbeam.gt.10) nbeam = 10 # MSR waiting for IMAS-3271
        resize!(dd.waves.coherent_wave,nbeam)
        for ibeam = 1:nbeam
            resize!(dd.waves.coherent_wave[ibeam].global_quantities, 1)
            dd.waves.coherent_wave[ibeam].identifier.antenna_name = dd.ec_launchers.beam[ibeam].name
            dd.waves.coherent_wave[ibeam].global_quantities[1].frequency = dd.ec_launchers.beam[ibeam].frequency.data[1]
            dd.waves.coherent_wave[ibeam].global_quantities[1].electrons.power_thermal = extrascal[ibeam, 2]
            dd.waves.coherent_wave[ibeam].global_quantities[1].power = extrascal[ibeam, 2]
            dd.waves.coherent_wave[ibeam].global_quantities[1].current_tor = extrascal[ibeam, 3]
            dd.waves.coherent_wave[ibeam].identifier.type.description = "TORBEAM"
            # if statement needed to filter empty sources for hcd2core_sources
            if (dd.ec_launchers.beam[ibeam].power_launched.data[1] > 0)
                dd.waves.coherent_wave[ibeam].identifier.type.name = "EC"
                dd.waves.coherent_wave[ibeam].identifier.type.index = 1
                dd.waves.coherent_wave[ibeam].wave_solver_type.index = 1 # BEAM/RAY TRACING
            end
            resize!(dd.waves.coherent_wave[ibeam].profiles_1d, 1)
            dd.waves.coherent_wave[ibeam].profiles_1d[1].time = @ddtime(dd.equilibrium.time)
            psi_beam = profout[ibeam, 1, 1:npnt].^2*(psiedge-psiax).+psiax
            rho_tor_beam = rho_tor_interpolator(psi_beam)
            dd.waves.coherent_wave[ibeam].profiles_1d[].grid.rho_tor_norm = rho_tor_beam
            dd.waves.coherent_wave[ibeam].profiles_1d[].grid.psi = psi_beam
            dd.waves.coherent_wave[ibeam].profiles_1d[].power_density = 1.e6*profout[ibeam, 2, 1:npnt]
            dd.waves.coherent_wave[ibeam].profiles_1d[].electrons.power_density_thermal = 1.e6*profout[ibeam, 2, 1:npnt]
            dd.waves.coherent_wave[ibeam].profiles_1d[].current_parallel_density = -1.e6*profout[ibeam, 3, 1:npnt]*iplasma/abs(iplasma)
            if ibeam == 1
                dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm = dd.waves.coherent_wave[1].profiles_1d[1].grid.rho_tor_norm
                dd.core_sources.source[source_index].profiles_1d[].grid.psi = dd.waves.coherent_wave[1].profiles_1d[1].grid.psi
                dd.core_sources.source[source_index].profiles_1d[].electrons.power_inside = zeros(Float64, length(dd.waves.coherent_wave[1].profiles_1d[].power_density))
                dd.core_sources.source[source_index].profiles_1d[].current_parallel_inside = zeros(Float64, length(dd.waves.coherent_wave[1].profiles_1d[].current_parallel_density))
            end
            dd.core_sources.source[source_index].profiles_1d[].electrons.power_inside += Interpolations.linear_interpolation(dd.waves.coherent_wave[ibeam].profiles_1d[1].grid.rho_tor_norm,
                                                dd.waves.coherent_wave[ibeam].profiles_1d[].power_density)(dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm)
            dd.core_sources.source[source_index].profiles_1d[].current_parallel_inside += Interpolations.linear_interpolation(dd.waves.coherent_wave[ibeam].profiles_1d[1].grid.rho_tor_norm,
                                                dd.waves.coherent_wave[ibeam].profiles_1d[].current_parallel_density)(dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm)
            
            
                                                # LOOP OVER RAYS
            resize!(dd.waves.coherent_wave[ibeam].beam_tracing, 1)
            dd.waves.coherent_wave[ibeam].beam_tracing[1].time = @ddtime(dd.equilibrium.time)
            resize!(dd.waves.coherent_wave[ibeam].beam_tracing[].beam,5) # Five beams/per gyrotron
            iend[] = min(npointsout[ibeam], ntraj)
            if (dd.ec_launchers.beam[ibeam].power_launched.data[1] > 0)
                for iray = 1:5 # 5 rays (to not mix with input beams)
                    r = 1.e-2*trajout[ibeam, 1+3*(iray-1), 1:iend[]]
                    z = 1.e-2*trajout[ibeam, 2+3*(iray-1), 1:iend[]]
                    phi = trajout[ibeam, 3+3*(iray-1), 1:iend[]] .+ dd.ec_launchers.beam[ibeam].launching_position.phi[1]
                    x = cos.(phi) .* r
                    y = sin.(phi) .* r
                    s = zeros(Float64, iend[])
                    s[2:iend[]] = sqrt.((x[2:iend[]] .- x[1:iend[]-1]).^2 .+ 
                                        (y[2:iend[]] .- y[1:iend[]-1]).^2 .+ 
                                        (z[2:iend[]] .- z[1:iend[]-1]).^2)
                    s = cumsum(s)
                    dd.waves.coherent_wave[ibeam].beam_tracing[].beam[iray].length = s
                    dd.waves.coherent_wave[ibeam].beam_tracing[].beam[iray].position.r = r
                    dd.waves.coherent_wave[ibeam].beam_tracing[].beam[iray].position.z = z
                    # Rotation of the output rays to fit the actual input toroidal angle
                    dd.waves.coherent_wave[ibeam].beam_tracing[].beam[iray].position.phi = phi
                end
            end

            

        end # LOOP OVER BEAMS (LAUNCHERS)
        dd.waves.code.output_flag = zeros(Int64,0) # NO ERROR

    end 

    const document = Dict()
    document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module TORBEAM
