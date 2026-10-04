module TORBEAM
using IMAS
import Libdl

Base.@kwdef struct TorbeamParams
    # switches
    npow::Int = 1             # Power absorption switch (1 = on, 0 = off)
    ncd::Int = 1              # Current drive calculation switch (1 = on, 0 = off)
    ncdroutine::Int = 2       # Current drive routine (0 = Curba, 1 = Lin-Liu, 2 = Lin-Liu + momentum conservation; Julia backend only: 3 = exact 2-D Lorentz-model solution, 4 = full linearized collision operator)
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

    # backend
    backend::Symbol = :fortran # :fortran calls libtorbeamB.so (needs TORBEAM_DIR); :julia is the pure-Julia implementation
end

# Array sizes from libtorbeam/src/libsrc/dimensions.f90
const MAXINT = 50
const MAXFLT = 50
const MAXRHR = 20
const NDAT = 100000
const NPNT = 5000
const MMAX = 450
const NMAX = MMAX
const MAXVOL = 100
const NTRAJ = 10000
const MAXDIM = 1 + MMAX + NMAX + 4 * MMAX * NMAX
const MAXLEN = 2 * MMAX + 2 * NMAX

"""
    BeamInputs

Everything the TORBEAM `beam` routine needs for one launcher: the integer and
float input vectors (`intinbeam`, `floatinbeam`), the equilibrium on the
rectangular R-Z grid (`eqdata`, `ni` x `nj`) and the ne/Te profiles (`prdata`,
`npsi` points). Units are TORBEAM's (cgs, MW, keV, 1e19 m^-3).
"""
struct BeamInputs
    intinbeam::Vector{Int32}
    floatinbeam::Vector{Float64}
    ni::Int
    nj::Int
    eqdata::Vector{Float64}
    npsi::Int
    prdata::Vector{Float64}
end

"""
    BeamOutputs

Raw results of one TORBEAM `beam` call, trimmed to their meaningful lengths:

  - `rhoresult`: scalars — `[1]` rho, `[2]` R, `[3]` Z at max absorption, `[13]` total driven current [kA],
    `[14]` total absorbed power [MW], `[20]` exit flag (0 = absorption, 1 = no plasma intersection,
    2 = crossed without absorption, 3 = cutoff, 4 = integrator failure, 5 = negative ne/Te, 6 = axis not found,
    7 = too many steps, 8 = cutoff at the vacuum-plasma boundary)
  - `t1data` (6*iend): R and Z of the central ray and of the upper/lower peripheral rays [cm]
  - `t1tdata` (6*iend): X and Y of the central ray and of the left/right peripheral rays [cm]
  - `t2data` (3*nprofv+9): flux-surface area/volume profile, then group velocity, widths and curvatures
  - `t2ndata` (3*NPNT): rho_pol, dP/dV [MW/m^3], j [MA/m^2]
  - `volprof` (2*MAXVOL): volume profile
"""
struct BeamOutputs
    rhoresult::Vector{Float64}
    iend::Int
    t1data::Vector{Float64}
    t1tdata::Vector{Float64}
    kend::Int
    t2data::Vector{Float64}
    t2ndata::Vector{Float64}
    icnt::Int
    ibgout::Int
    volprof::Vector{Float64}
end

include("model.jl")
include("dispersion.jl")
include("absorption.jl")
include("beam_tracing.jl")
include("spitzer.jl")
include("current_drive.jl")
include("deposition.jl")
include("julia_backend.jl")

"""
    equilibrium_inputs(dd::IMAS.dd)

Assemble the launcher-independent TORBEAM inputs from the current equilibrium
time slice and core profiles: the `eqdata` and `prdata` vectors plus the
grid/profile sizes and the sign of psi at the axis.
"""
function equilibrium_inputs(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
    cp1d = dd.core_profiles.profiles_1d[]

    eqdata = zeros(Float64, MAXDIM)
    prdata = zeros(Float64, MAXLEN)

    # data in IMAS format
    Rarr = eqt2d.grid.dim1
    Zarr = eqt2d.grid.dim2
    ni = length(Rarr)
    nj = length(Zarr)
    br = eqt2d.b_field_r
    bt = eqt2d.b_field_tor
    bz = eqt2d.b_field_z

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

    # DETERMINE WHETHER PSI FLUX IS MAXIMUM (1) OR MINIMUM (-1) AT THE MAGNETIC AXIS
    sgnm = psiedge > psiax ? 1.0 : -1.0

    return (; eqdata, ni, nj, prdata, npsi, sgnm, psiedge, psiax)
end

"""
    beam_inputs(dd::IMAS.dd, ibeam::Int, torbeam_params::TorbeamParams, eq)

Assemble the `BeamInputs` for launcher `ibeam`, given the launcher-independent
part `eq` from [`equilibrium_inputs`](@ref).
"""
function beam_inputs(dd::IMAS.dd, ibeam::Int, torbeam_params::TorbeamParams, eq)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    beam = dd.ec_launchers.beam[ibeam]
    ps_beam = dd.pulse_schedule.ec.beam[ibeam]
    power_launched = @ddtime(ps_beam.power_launched.reference)

    intinbeam = zeros(Int32, MAXINT)
    floatinbeam = zeros(Float64, MAXFLT)

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
    intinbeam[7] = torbeam_params.ncdroutine > 2 ? torbeam_params.ncdroutine - 2 : torbeam_params.ncdroutine   # Julia-only 3/4 → nearest Fortran model
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
    intinbeam[MAXINT] = torbeam_params.nrel

    #floatinbeam(17:18): obsolete --> not filled)
    #floatinbeam(6:13):  analytic --> not filled)
    #floatinbeam(26:32): analytic --> not filled)
    floatinbeam[1] = @ddtime(beam.frequency.data)  # (xf)
    # floatinbeam[2] = rad2deg(-@ddtime(beam.steering_angle_tor))
    # floatinbeam[3] = rad2deg(@ddtime(beam.steering_angle_pol))
    # TODO fix when OMAS is updated
    steering_angle_tor = -asin(cos(@ddtime(beam.steering_angle_pol)) * sin(@ddtime(beam.steering_angle_tor)))
    steering_angle_pol = atan(tan(@ddtime(beam.steering_angle_pol)), cos(@ddtime(beam.steering_angle_tor)))
    alpha = steering_angle_pol
    beta = -steering_angle_tor
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
    floatinbeam[34] = eq.sgnm  # (deduced from psi_ed-psi_ax)
    floatinbeam[35] = cp1d.zeff[1]  # (xzeff)
    floatinbeam[36] = torbeam_params.rhostop  # keep
    floatinbeam[37] = torbeam_params.xzsrch  # keep

    return BeamInputs(intinbeam, floatinbeam, eq.ni, eq.nj, eq.eqdata, eq.npsi, eq.prdata)
end

"""
    run_beam(inputs::BeamInputs, torbeam_params::TorbeamParams; cache=nothing)

Run TORBEAM for one launcher with the backend selected in `torbeam_params`.
`cache` (a `BackendCache`, see `backend_cache`) is used by the Julia backend to
share the per-equilibrium objects between the beams of one run.
"""
function run_beam(inputs::BeamInputs, torbeam_params::TorbeamParams; cache=nothing)
    if torbeam_params.backend == :fortran
        return fortran_beam(inputs, torbeam_params)
    elseif torbeam_params.backend == :julia
        return julia_beam(inputs, torbeam_params; cache)
    else
        error("TORBEAM backend `$(torbeam_params.backend)` not implemented (available: :fortran, :julia)")
    end
end

"""
    fortran_library()

Path of `libtorbeamB.so` as resolved from `TORBEAM_DIR`
"""
fortran_library() = get(ENV, "TORBEAM_DIR", "") * "/../lib/libtorbeamB.so"

"""
    fortran_available()

Whether the TORBEAM Fortran library can be found (`TORBEAM_DIR` set and the `.so` present)
"""
fortran_available() = haskey(ENV, "TORBEAM_DIR") && isfile(fortran_library())

"""
    fortran_beam(inputs::BeamInputs, torbeam_params::TorbeamParams)

Call the `beam` routine of `libtorbeamB.so` and return its trimmed outputs
"""
function fortran_beam(inputs::BeamInputs, torbeam_params::TorbeamParams)
    nprofv = torbeam_params.nprofv

    # Define output scalars
    iend = Ref{Int32}(0)
    kend = Ref{Int32}(0)
    icnt = Ref{Int32}(0)
    ibgout = Ref{Int32}(0)

    # Allocate output arrays
    rhoresult = zeros(Float64, MAXRHR)
    t1data = zeros(Float64, 6 * NDAT)
    t1tdata = zeros(Float64, 6 * NDAT)
    t2data = zeros(Float64, 5 * NDAT)
    t2ndata = zeros(Float64, 3 * NPNT)
    volprof = zeros(Float64, 2 * MAXVOL)

    function invoke_ccall()
        # The library path is only known at run time (TORBEAM_DIR), so resolve the
        # symbol through Libdl: `ccall((:sym, <non-constant expr>), ...)` is rejected
        # at lowering time since Julia 1.13.
        beam_ptr = Libdl.dlsym(Libdl.dlopen(fortran_library()), :beam_)   # Name in the shared library (append `_`)
        return ccall(
            beam_ptr,
            Cvoid,                             # Return type
            (Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, # Inputs
                Ptr{Float64}, Ref{Cint}, Ptr{Float64}, Ptr{Float64}, Ref{Cint},
                Ptr{Float64}, Ptr{Float64}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{Float64}), # Argument types
            inputs.intinbeam, inputs.floatinbeam, inputs.ni, inputs.nj, inputs.eqdata, inputs.npsi, inputs.npsi, inputs.prdata,
            rhoresult, iend, t1data, t1tdata, kend,
            t2data, t2ndata, icnt, ibgout, nprofv, volprof
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

    return BeamOutputs(
        rhoresult,
        iend[],
        t1data[1:6*iend[]],
        t1tdata[1:6*iend[]],
        kend[],
        t2data[1:3*nprofv+9],
        t2ndata,
        icnt[],
        ibgout[],
        volprof)
end

"""
    ray_trajectories(out::BeamOutputs)

Central ray and 4 peripheral rays as a `(15, npoints)` matrix of
`(r [cm], z [cm], phi [rad])` triplets, one triplet per ray.
"""
function ray_trajectories(out::BeamOutputs)
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
    iend = out.iend
    t1data = out.t1data
    t1tdata = out.t1tdata

    # TO PREVENT GOING BEYOND THE PRE-DEFINED TRAJOUT ARRAY
    npoints = min(iend, NTRAJ)
    trajout = zeros(Float64, 15, npoints)

    # TRAJECTORY OF CENTRAL RAY AND 4 PERIPHERAL "RAYS"
    for lfd in 1:npoints
        # 1st ray
        trajout[1, lfd] = t1data[lfd]                    # r
        trajout[2, lfd] = t1data[iend+lfd]               # z
        trajout[3, lfd] = atan(t1tdata[iend+lfd], t1tdata[lfd]) # phi
        # 2nd ray
        trajout[4, lfd] = t1data[2*iend+lfd]
        trajout[5, lfd] = t1data[3*iend+lfd]
        trajout[6, lfd] = atan(t1tdata[iend+lfd], t1tdata[lfd])
        # 3rd ray
        trajout[7, lfd] = t1data[4*iend+lfd]
        trajout[8, lfd] = t1data[5*iend+lfd]
        trajout[9, lfd] = atan(t1tdata[iend+lfd], t1tdata[lfd])
        # 4th ray
        trajout[10, lfd] = sqrt(t1tdata[2*iend+lfd]^2 + t1tdata[3*iend+lfd]^2)
        trajout[11, lfd] = t1data[iend+lfd]
        trajout[12, lfd] = atan(t1tdata[3*iend+lfd], t1tdata[2*iend+lfd])
        # 5th ray
        trajout[13, lfd] = sqrt(t1tdata[4*iend+lfd]^2 + t1tdata[5*iend+lfd]^2)
        trajout[14, lfd] = t1data[iend+lfd]
        trajout[15, lfd] = atan(t1tdata[5*iend+lfd], t1tdata[4*iend+lfd])
    end

    return trajout
end

"""
    run_torbeam(dd::IMAS.dd, torbeam_params::TorbeamParams)

Run TORBEAM for all launchers with non-zero power at the current time and store
the results in the `waves` and `core_sources` IDSs
"""
function run_torbeam(dd::IMAS.dd, torbeam_params::TorbeamParams)
    nbeam = length(dd.ec_launchers.beam)
    if nbeam < 1
        return nbeam
    end

    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # Interpolator for psi -> rho-tor need that later
    rho_tor_norm_interpolator = IMAS.interp1d(eqt1d.psi, eqt1d.rho_tor_norm)

    eq = equilibrium_inputs(dd)
    psiedge = eq.psiedge
    psiax = eq.psiax

    outputs = Vector{Union{Nothing,BeamOutputs}}(nothing, nbeam)
    cache = nothing   # per-equilibrium objects of the Julia backend, built on the first active beam
    # LOOP OVER BEAMS OF THE EC_LAUNCHERS IDS
    for ibeam in 1:nbeam
        ps_beam = dd.pulse_schedule.ec.beam[ibeam]
        power_launched = @ddtime(ps_beam.power_launched.reference)

        # ONLY DEAL WITH ACTIVE BEAMS
        if power_launched > 0
            inputs = beam_inputs(dd, ibeam, torbeam_params, eq)
            if torbeam_params.backend == :julia && cache === nothing
                cache = backend_cache(inputs, torbeam_params)
            end
            @debug("------------------------------------------------------------")
            @debug("Input power beam: ", ibeam, " ", power_launched * 1.e-6, " MW")
            outputs[ibeam] = run_beam(inputs, torbeam_params; cache)
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
        out = outputs[ibeam]

        wv = dd.waves.coherent_wave[ibeam]
        wv.identifier.antenna_name = beam.name
        wv.identifier.type.description = "TORBEAM"
        wv.identifier.type.name = "EC"
        wv.identifier.type.index = 1
        wv.wave_solver_type.index = 1 # BEAM/RAY TRACING

        wvg = resize!(wv.global_quantities) # global_time
        wvg.frequency = @ddtime(beam.frequency.data)
        # rhoresult(13) [MW] and rhoresult(12) [kA] in TORBEAM's 0-based numbering
        wvg.electrons.power_thermal = out === nothing ? 0.0 : 1.e6 * out.rhoresult[14]
        wvg.power = wvg.electrons.power_thermal
        wvg.current_tor = out === nothing ? 0.0 : 1.e3 * out.rhoresult[13]

        # 1D PROFILES OF RHO, DP/DV, J (CHECKED OK, DIM = NPNT = 5000)
        wv1d = resize!(wv.profiles_1d) # global_time
        wv.profiles_1d[1].time = @ddtime(dd.equilibrium.time)
        rho_pol = out === nothing ? zeros(NPNT) : out.t2ndata[1:NPNT]
        dPdV = out === nothing ? zeros(NPNT) : out.t2ndata[NPNT+1:2*NPNT]
        j = out === nothing ? zeros(NPNT) : out.t2ndata[2*NPNT+1:3*NPNT]
        psi_beam = rho_pol .^ 2 * (psiedge - psiax) .+ psiax
        rho_tor_norm_beam = rho_tor_norm_interpolator.(psi_beam)
        wv1d.grid.rho_tor_norm = rho_tor_norm_beam
        wv1d.grid.psi = psi_beam
        wv1d.power_density = 1.e6 * dPdV
        wv1d.electrons.power_density_thermal = 1.e6 * dPdV
        wv1d.current_parallel_density = -1.e6 * j * sign(eqt.global_quantities.ip)

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
        if out !== nothing
            trajout = ray_trajectories(out)
            npoints = size(trajout, 2)
            for iray in 1:torbeam_params.n_ray
                r = 1.e-2 * trajout[1+3*(iray-1), :]
                z = 1.e-2 * trajout[2+3*(iray-1), :]
                # FIX after OMAS ec_launchers correction
                phi_launch = -beam.launching_position.phi[1] - pi / 2.0
                phi = trajout[3+3*(iray-1), :] .+ phi_launch
                x = cos.(phi) .* r
                y = sin.(phi) .* r
                s = zeros(Float64, npoints)
                s[2:npoints] = sqrt.((x[2:npoints] .- x[1:npoints-1]) .^ 2 .+
                                     (y[2:npoints] .- y[1:npoints-1]) .^ 2 .+
                                     (z[2:npoints] .- z[1:npoints-1]) .^ 2)
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

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
