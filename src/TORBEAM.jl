module TORBEAM
using IMAS

Base.@kwdef struct TorbeamParams
    # switches
    npow::Int = 1             # Power absorption switch (1 = on, 0 = off)
    ncd::Int = 1              # Current drive calculation switch (1 = on, 0 = off)
    ncdroutine::Int = 2       # Current drive routine selection (0 = Curba, 1 = Lin-Liu, 2 = Lin-Liu + momentum conservation)
    nprofv::Int = 50          # Number of radial points for volume profile calculation
    noout::Int = 0            # Screen output switch (0 = output enabled, 1 = output disabled)
    nrela::Int = 1            # Relativity consideration in absorption (0 = weakly, 1 = fully relativistic)
    nmaxh::Int = 3            # Number of harmonics to consider (1 to 5)
    nabsroutine::Int = 1      # Absorption routine selection (0 = Westerhof, 1 = Farina)
    nastra::Int = 0           # Definition of driven current density (0 = Lin-Liu, 1 = ASTRA, 2 = JINTRAC)
    nprofcalc::Int = 1        # Deposition profile calculation method (0 = standard, 1 = Maj method)
    ncdharm::Int = 1          # Harmonic consideration in current drive efficiency (0 = lowest harmonic only, 1 = includes next harmonic)
    nrel::Int = 0             # Relativistic mass correction for reflectometry (1 = enabled, 0 = disabled)
    n_ray::Int = 5            # Number of rays used in beam tracing
    verbose::Bool = false     # Verbose mode for output debugging (true = enabled, false = disabled)

    # Float parameters
    xrtol::Float64 = 1e-07    # Required relative error tolerance
    xatol::Float64 = 1e-07    # Required absolute error tolerance
    xstep::Float64 = 2.0      # Integration step in vacuum (cm)
    rhostop::Float64 = 0.96   # Maximum value of the flux coordinate (rho) before stopping
    xzsrch::Float64 = 0.0     # Vertical position for searching the magnetic axis (default 0 cm)
end

function run_torbeam(dd::IMAS.dd, torbeam_params::TorbeamParams)
    # TORBEAM subroutine in IMAS
    #----------------------------

    # Define parameters
    # from libtorbeam/src/libsrc/dimensions.f90 -> Could be moved to a torbeam.yaml (?)
    # Define sizes of various arrays. Need to allocate them here to then pass to TORBEAM

    nbeam = length(dd.ec_launchers.beam)
    if nbeam < 1
        return nbeam
    end

    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
    cp1d = dd.core_profiles.profiles_1d[]

    maxint = 50
    maxflt = 50
    maxrhr = 20
    ndat = 100000
    npnt = 5000
    mmax = 450
    maxvol = 100
    ntraj = 10000
    nmax = mmax
    maxdim = 1 + mmax + nmax + 4 * mmax * nmax
    maxlen = 2 * mmax + 2 * nmax

    # Allocate input arrays
    rhoresult = zeros(Float64, maxrhr)
    intinbeam = zeros(Int32, maxint)
    floatinbeam = zeros(Float64, maxflt)
    eqdata = zeros(Float64, maxdim)
    prdata = zeros(Float64, maxlen)
    volprof = zeros(Float64, 2 * maxvol)

    # Define output scalars
    iend = Ref{Int32}(0)
    kend = Ref{Int32}(0)
    icnt = Ref{Int32}(0)
    ibgout = Ref{Int32}(0)

    # Allocate output arrays
    t1data = zeros(Float64, 6 * ndat)
    t1tdata = zeros(Float64, 6 * ndat)
    t2data = zeros(Float64, 5 * ndat)
    t2ndata = zeros(Float64, 3 * npnt)
    trajout = zeros(Float64, (nbeam, 15, ntraj))

    # Allocate arrays
    Rarr = eqt2d.grid.dim1
    Zarr = eqt2d.grid.dim2
    ni = length(Rarr)
    nj = length(Zarr)
    br = eqt2d.b_field_r
    bt = eqt2d.b_field_tor
    bz = eqt2d.b_field_z

    # Interpolator for psi -> rho-tor need that later
    rho_tor_norm_interpolator = IMAS.interp1d(eq1d.psi, eq1d.rho_tor_norm)

    psiedge = eqt.global_quantities.psi_boundary
    psiax = eqt.global_quantities.psi_axis

    # Initialize vector 'eqdata' as TORBEAM input (topfile):
    # Psi, 1D, for normalization
    npsi = length(cp1d.grid.psi)
    psi = cp1d.grid.psi
    # To ensure consistency between the 1D and 2D psi profiles: take both from the equilibrium IDS

    # FILL TORBEAM INTERNAL EQUILIBRIUM DATA
    eqdata[1] = psiedge
    for i in 1:ni
        eqdata[i+1] = Rarr[i]
    end
    for j in 1:nj
        eqdata[j+ni+1] = Zarr[j]
    end
    k = 0
    for j in 1:nj
        for i in 1:ni
            k = k + 1
            eqdata[k+ni+nj+1] = br[i, j]
        end
    end
    k = 0
    for j in 1:nj
        for i in 1:ni
            k = k + 1
            eqdata[k+ni+nj+ni*nj+1] = bt[i, j]
        end
    end
    k = 0
    for j in 1:nj
        for i in 1:ni
            k = k + 1
            eqdata[k+ni+nj+2*ni*nj+1] = bz[i, j]
        end
    end
    k = 0
    for j in 1:nj
        for i in 1:ni
            k = k + 1
            eqdata[k+ni+nj+3*ni*nj+1] = eqt2d.psi[i, j]
        end
    end

    # FILL TORBEAM INTERNAL PROFILE DATA
    #... Initialize vector 'prdata' as TORBEAM input (ne.dat & Te.dat):
    # Psi and profiles
    for i in 1:npsi
        prdata[i] = sqrt((psi[i] - psi[1]) / (psi[npsi] - psi[1]))
        prdata[i+npsi] = cp1d.electrons.density[i] * 1.0e-19
    end
    for i in 1:npsi
        prdata[i+2*npsi] = sqrt((psi[i] - psi[1]) / (psi[npsi] - psi[1]))
        prdata[i+2*npsi+npsi] = cp1d.electrons.temperature[i] * 1.0e-3
    end

    # Initialize antenna data as TORBEAM input:
    # DETERMINE WHETHER PSI FLUX IS MAXIMUM (1) OR MINIMUM (-1) AT THE MAGNETIC AXIS
    if (psiedge > psiax)
        sgnm = 1.0
    else
        sgnm = -1.0
    end
    npointsout = zeros(Int64, nbeam)
    extrascal = zeros(Float64, (nbeam, 5))
    profout = zeros(Float64, (nbeam, 3, npnt))
    # LOOP OVER BEAMS OF THE EC_LAUNCHERS IDS
    for ibeam in 1:nbeam
        beam = dd.ec_launchers.beam[ibeam]
        ps_beam = dd.pulse_schedule.ec.beam[ibeam]
        power_launched = @ddtime(ps_beam.power_launched.reference)

        # ONLY DEAL WITH ACTIVE BEAMS
        if power_launched > 0

            # IT LOOKS LIKE TORBEAM NEEDS PHI = 0, OTHERWISE IT DOES NOT TREAT THE BEAM PROPERLY
            # BUT WE WILL RESTORE THE ACTUAL PHI ANGLE AFTER THE RAY-TRACING, SO WE DON'T
            # PUT ec_launchers%BEAM(IBEAM)%LAUNCHING_POSITION%PHI TO 0 ANYMORE
            # (WE ARTIFICIALLY PUT PHI=0 IN floatinbeam(3) AND FLOTINBEAM(4) INSTEAD

            #intinbeam
            intinbeam[1] = 2  # tbr
            intinbeam[2] = 2  # tbr
            intinbeam[3] = beam.mode  # (nmod)
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
            floatinbeam[1] = @ddtime(beam.frequency.data)  # (xf)
            alpha = @ddtime(dd.ec_launchers.beam[ibeam].steering_angle_pol)
            beta = -@ddtime(dd.ec_launchers.beam[ibeam].steering_angle_tor)
            floatinbeam[2] = rad2deg(atan(tan(beta), cos(alpha)))
            floatinbeam[3] = rad2deg(asin(sin(alpha) * cos(beta)))
            floatinbeam[4] = 1.e2 * beam.launching_position.r[1] * cos(0)  # (xxb)
            floatinbeam[5] = 1.e2 * beam.launching_position.r[1] * sin(0)  # (xyb)
            floatinbeam[6] = 1.e2 * beam.launching_position.z[1]  # (xzb)

            floatinbeam[15] = torbeam_params.xrtol  # keep
            floatinbeam[16] = torbeam_params.xatol  # keep
            floatinbeam[17] = torbeam_params.xstep  # keep
            floatinbeam[20] = -1.e2 / (beam.phase.curvature[1, 1])  # (xryyb)
            floatinbeam[21] = -1.e2 / (beam.phase.curvature[2, 1])  # (xrzzb)
            if (cos(@ddtime(beam.spot.angle))^2 > 0.5)
                floatinbeam[22] = beam.spot.size[1, 1] * 1.e2  # (xwyyb)
                floatinbeam[23] = beam.spot.size[2, 1] * 1.e2  # (xwzzb)
            else
                floatinbeam[22] = beam.spot.size[2, 1] * 1.e2  # (xwzzb)
                floatinbeam[23] = beam.spot.size[1, 1] * 1.e2  # (xwyyb)
            end
            floatinbeam[24] = power_launched * 1.e-6  # (xpw0)
            floatinbeam[25] = eqt.boundary.geometric_axis.r * 1e2  # (xrmaj)
            floatinbeam[26] = eqt.boundary.minor_radius * 1e2  # (xrmin)
            floatinbeam[27] = eqt.global_quantities.vacuum_toroidal_field.b0
            floatinbeam[34] = sgnm  # (deduced from psi_ed-psi_ax)
            floatinbeam[35] = cp1d.zeff[1]  # (xzeff)
            floatinbeam[36] = torbeam_params.rhostop  # keep
            floatinbeam[37] = torbeam_params.xzsrch  # keep

            @debug("------------------------------------------------------------")
            @debug("Input power beam: ", ibeam, " ", power_launched * 1.e-6, " MW")
            # CALL TORBEAM
            function invoke_ccall()
                return ccall(
                    (:beam_, get(ENV, "TORBEAM_DIR", "") * "/../lib/libtorbeamB.so"),    # Name in the shared library (append `_`)
                    Cvoid,                             # Return type
                    (Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, # Inputs
                        Ptr{Float64}, Ref{Cint}, Ptr{Float64}, Ptr{Float64}, Ref{Cint},
                        Ptr{Float64}, Ptr{Float64}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{Float64}), # Argument types
                    intinbeam, floatinbeam, ni, nj, eqdata, npsi, npsi, prdata,
                    rhoresult, iend, t1data, t1tdata, kend,
                    t2data, t2ndata, icnt, ibgout, torbeam_params.nprofv, volprof
                )
            end

            if torbeam_params.verbose
                invoke_ccall()
            else
                redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        return invoke_ccall()
                    end
                end
            end

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
            extrascal[ibeam, 1] = power_launched * 1.e-6
            extrascal[ibeam, 2] = 1.e6 * rhoresult[13]
            extrascal[ibeam, 3] = 1.e6 * rhoresult[12]

            # 1D PROFILES OF RHO, DP/DV, J (CHECKED OK, DIM = NPNT = 5000)
            for kp in 0:2
                for lfd in 1:npnt
                    profout[ibeam, kp+1, lfd] = t2ndata[kp*npnt+lfd]
                end
            end

            # TO PREVENT GOING BEYOND THE PRE-DEFINED TRAJOUT ARRAY
            if (iend[] > ntraj)
                iend[] = ntraj
            end

            # TRAJECTORY OF CENTRAL RAY AND 4 PERIPHERAL "RAYS"
            for lfd in 1:iend[]
                # 1st ray
                trajout[ibeam, 1, lfd] = t1data[lfd]                    # r
                trajout[ibeam, 2, lfd] = t1data[iend[]+lfd]               # z
                trajout[ibeam, 3, lfd] = atan(t1tdata[iend[]+lfd], t1tdata[lfd]) # phi
                # 2nd ray
                trajout[ibeam, 4, lfd] = t1data[2*iend[]+lfd]
                trajout[ibeam, 5, lfd] = t1data[3*iend[]+lfd]
                trajout[ibeam, 6, lfd] = atan(t1tdata[iend[]+lfd], t1tdata[lfd])
                # 3rd ray
                trajout[ibeam, 7, lfd] = t1data[4*iend[]+lfd]
                trajout[ibeam, 8, lfd] = t1data[5*iend[]+lfd]
                trajout[ibeam, 9, lfd] = atan(t1tdata[iend[]+lfd], t1tdata[lfd])
                # 4th ray
                trajout[ibeam, 10, lfd] = sqrt(t1tdata[2*iend[]+lfd]^2 + t1tdata[3*iend[]+lfd]^2)
                trajout[ibeam, 11, lfd] = t1data[iend[]+lfd]
                trajout[ibeam, 12, lfd] = atan(t1tdata[3*iend[]+lfd], t1tdata[2*iend[]+lfd])
                # 5th ray
                trajout[ibeam, 13, lfd] = sqrt(t1tdata[4*iend[]+lfd]^2 + t1tdata[5*iend[]+lfd]^2)
                trajout[ibeam, 14, lfd] = t1data[iend[]+lfd]
                trajout[ibeam, 15, lfd] = atan(t1tdata[5*iend[]+lfd], t1tdata[4*iend[]+lfd])
            end

        end # TEST BEAM_POWER > 0

    end # LOOP OVER BEAMS OF EC_LAUNCHERS IDS

    # ----------------------------
    # SAVE RESULTS INTO WAVES IDS
    # ----------------------------

    # LOOP OVER BEAMS (LAUNCHERS)
    #if(nbeam.gt.10) nbeam = 10 # MSR waiting for IMAS-3271
    resize!(dd.waves.coherent_wave, nbeam)
    for ibeam in 1:nbeam
        beam = dd.ec_launchers.beam[ibeam]
        ps_beam = dd.pulse_schedule.ec.beam[ibeam]
        power_launched = @ddtime(ps_beam.power_launched.reference)

        wv = dd.waves.coherent_wave[ibeam]
        wv.identifier.antenna_name = beam.name
        wv.identifier.type.description = "TORBEAM"
        wv.identifier.type.name = "EC"
        wv.identifier.type.index = 1
        wv.wave_solver_type.index = 1 # BEAM/RAY TRACING

        wvg = resize!(wv.global_quantities) # global_time
        wvg.frequency = @ddtime(beam.frequency.data)
        wvg.electrons.power_thermal = extrascal[ibeam, 2]
        wvg.power = extrascal[ibeam, 2]
        wvg.current_tor = extrascal[ibeam, 3]

        wv1d = resize!(wv.profiles_1d) # global_time
        psi_beam = profout[ibeam, 1, 1:npnt] .^ 2 * (psiedge - psiax) .+ psiax
        rho_tor_norm_beam = rho_tor_norm_interpolator.(psi_beam)
        wv1d.grid.rho_tor_norm = rho_tor_norm_beam
        wv1d.grid.psi = psi_beam
        wv1d.power_density = 1.e6 * profout[ibeam, 2, 1:npnt]
        wv1d.electrons.power_density_thermal = 1.e6 * profout[ibeam, 2, 1:npnt]
        #TODO: The expression below might still need to be revised
        # wv1d.current_parallel_density = -1.e6 * profout[ibeam, 3, 1:npnt] * sign(eqt.global_quantities.ip)
        wv1d.current_parallel_density = 1.e6 * profout[ibeam, 3, 1:npnt]
        
        source = resize!(dd.core_sources.source, :ec, "identifier.name" => beam.name; wipe=false)
        IMAS.new_source(
            source,
            source.identifier.index,
            beam.name,
            wv1d.grid.rho_tor_norm,
            wv1d.grid.volume,
            wv1d.grid.area;
            electrons_energy=wv1d.power_density,
            j_parallel=wv1d.current_parallel_density)

        # LOOP OVER RAYS
        wvb = resize!(wv.beam_tracing) # global_time
        resize!(wvb.beam, torbeam_params.n_ray) # Five beams/per gyrotron
        iend[] = min(npointsout[ibeam], ntraj)
        if power_launched > 0.0
            for iray in 1:torbeam_params.n_ray
                r = 1.e-2 * trajout[ibeam, 1+3*(iray-1), 1:iend[]]
                z = 1.e-2 * trajout[ibeam, 2+3*(iray-1), 1:iend[]]
                phi_launch = beam.launching_position.phi[1]
                phi = trajout[ibeam, 3+3*(iray-1), 1:iend[]] .+ phi_launch
                x = cos.(phi) .* r
                y = sin.(phi) .* r
                s = zeros(Float64, iend[])
                s[2:iend[]] = sqrt.((x[2:iend[]] .- x[1:iend[]-1]) .^ 2 .+
                                    (y[2:iend[]] .- y[1:iend[]-1]) .^ 2 .+
                                    (z[2:iend[]] .- z[1:iend[]-1]) .^ 2)
                s = cumsum(s)
                wvb.beam[iray].length = s
                wvb.beam[iray].position.r = r
                wvb.beam[iray].position.z = z
                # Rotation of the output rays to fit the actual input toroidal angle
                wvb.beam[iray].position.phi = phi
            end
        end

    end # LOOP OVER BEAMS (LAUNCHERS)

    @ddtime(dd.waves.code.output_flag = 0) # NO ERROR

    return nbeam
end

function overview_plot end



const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
