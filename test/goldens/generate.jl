# Generate the regression ("golden") data for TORBEAM.jl from the Fortran library.
#
# For each case this script
#   1. builds a `dd` with FUSE (equilibrium, core profiles, EC launcher),
#   2. replaces the EC launcher with an explicit, physically sensible set of
#      beams (frequency, mode, launch point, steering angles, focusing, power) —
#      FUSE's automatic launcher placement is not usable here: it puts the
#      launcher in the wall corner and, for DIII-D, picks a frequency in cutoff,
#   3. stores the trimmed `dd` in test/data/<case>.json,
#   4. runs the Fortran `beam` routine on every beam and stores its raw outputs
#      (BeamOutputs) in test/goldens/<case>.json.
#
# Needs FUSE and the Fortran library, so it is NOT part of `Pkg.test`. On omega:
#
#   module load torbeam/gcc11.x
#   julia --project=@torbeam-dev test/goldens/generate.jl      # env with FUSE, JSON and this TORBEAM dev'ed
#
# The tests (test/runtests.jl) then only need the JSON files, and skip the
# Fortran comparison when the library is absent.

using FUSE
using IMAS
using TORBEAM
using JSON
using Dates

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const GOLDENS_DIR = @__DIR__

@assert TORBEAM.fortran_available() "TORBEAM Fortran library not found: `module load torbeam` first (TORBEAM_DIR=$(get(ENV, "TORBEAM_DIR", "")))"

const FOCUS = -1 / 1.5  # phase curvature [1/m] giving TORBEAM a +150 cm wavefront radius

"""
Each case: how to build the `dd`, the launcher (frequency [Hz], mode, launch R/Z [m])
and the beam variants (`pol`/`tor` steering angles [deg], phase `curvature` [1/m],
`power` fraction of the case's launched power, optional `mode` override).
"""
const CASES = Dict(
    "D3D" => (
        build=() -> FUSE.init(FUSE.case_parameters(:D3D, :default)...),
        frequency=110e9, mode=-1, launch=(2.40, 0.67),     # DIII-D 110 GHz gyrotrons, X2 at 1.7 T
        variants=[
            (name="pol20", pol=20.0, tor=0.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol20_focus_flipped", pol=20.0, tor=0.0, curvature=-FOCUS, power=1.0, mode=nothing),
            (name="pol10_tor15", pol=10.0, tor=15.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol30_tor-15", pol=30.0, tor=-15.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol30_tor15_half_power", pol=30.0, tor=15.0, curvature=FOCUS, power=0.5, mode=nothing),
            (name="pol30_tor15_Omode", pol=30.0, tor=15.0, curvature=FOCUS, power=1.0, mode=1),
        ]),
    "ITER" => (
        build=() -> FUSE.init(FUSE.case_parameters(:ITER; init_from=:ods)...),
        frequency=138.134715470e9, mode=1, launch=(8.39, 4.72),  # O1 upper launcher
        variants=[
            (name="pol63", pol=63.0, tor=0.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol63_focus_flipped", pol=63.0, tor=0.0, curvature=-FOCUS, power=1.0, mode=nothing),
            (name="pol60_tor20", pol=60.0, tor=20.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol66_tor-20", pol=66.0, tor=-20.0, curvature=FOCUS, power=1.0, mode=nothing),
            (name="pol63_tor20_half_power", pol=63.0, tor=20.0, curvature=FOCUS, power=0.5, mode=nothing),
            (name="pol66", pol=66.0, tor=0.0, curvature=FOCUS, power=1.0, mode=nothing),
        ]),
)

"""
    make_beams!(dd, spec)

Replace `dd.ec_launchers.beam` (and the matching pulse_schedule beams) with one
beam per variant in `spec`, all derived from the first launcher.
"""
function make_beams!(dd::IMAS.dd, spec)
    # read the baseline while the IDSs are still attached to dd (time coordinates
    # such as pulse_schedule.ec.time live on the parents)
    power0 = @ddtime(dd.pulse_schedule.ec.beam[1].power_launched.reference)
    beam0 = deepcopy(dd.ec_launchers.beam[1])
    ps0 = deepcopy(dd.pulse_schedule.ec.beam[1])
    n = length(spec.variants)
    resize!(dd.ec_launchers.beam, n; wipe=true)
    resize!(dd.pulse_schedule.ec.beam, n; wipe=true)
    for (k, v) in enumerate(spec.variants)
        # setindex! re-attaches the copies to dd
        beam = dd.ec_launchers.beam[k] = deepcopy(beam0)
        ps = dd.pulse_schedule.ec.beam[k] = deepcopy(ps0)
        beam.name = ps.name = "ec_$(v.name)"
        @ddtime(beam.frequency.data = spec.frequency)
        beam.mode = v.mode === nothing ? spec.mode : v.mode
        @ddtime(beam.launching_position.r = spec.launch[1])
        @ddtime(beam.launching_position.z = spec.launch[2])
        @ddtime(beam.launching_position.phi = 0.0)
        @ddtime(beam.steering_angle_pol = deg2rad(v.pol))
        @ddtime(beam.steering_angle_tor = deg2rad(v.tor))
        @ddtime(beam.phase.curvature = [v.curvature, v.curvature])
        @ddtime(ps.power_launched.reference = power0 * v.power)
    end
    return dd
end

"""
    trimmed(dd)

Only the IDSs TORBEAM reads, with expressions frozen, so the JSON is self-contained
"""
function trimmed(dd::IMAS.dd)
    dd2 = IMAS.dd()
    dd2.global_time = dd.global_time
    dd2.equilibrium = deepcopy(dd.equilibrium)
    dd2.core_profiles = deepcopy(dd.core_profiles)
    dd2.ec_launchers = deepcopy(dd.ec_launchers)
    dd2.pulse_schedule.ec = deepcopy(dd.pulse_schedule.ec)
    dd2.wall = deepcopy(dd.wall)
    return dd2
end

# 10 significant digits are plenty for the 1e-6 test tolerances and keep the files small
round10(x::Float64) = round(x; sigdigits=10)
round10(x::AbstractArray) = round10.(x)
round10(x::Dict) = Dict(k => round10(v) for (k, v) in x)
round10(x::AbstractVector{Any}) = [round10(v) for v in x]
round10(x) = x

params = TORBEAM.TorbeamParams()
params_dict = Dict(string(f) => getfield(params, f) for f in fieldnames(TORBEAM.TorbeamParams))

for (case, spec) in CASES
    println("### $case")
    dd = spec.build()
    make_beams!(dd, spec)
    dd = trimmed(dd)
    IMAS.imas2json(dd, joinpath(DATA_DIR, "$case.json"); freeze=true)

    eq = TORBEAM.equilibrium_inputs(dd)
    beams = []
    for ibeam in eachindex(dd.ec_launchers.beam)
        inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
        out = TORBEAM.fortran_beam(inputs, params)
        println("  beam $ibeam $(dd.ec_launchers.beam[ibeam].name): P_abs = $(round(out.rhoresult[14]; digits=4)) MW, I_cd = $(round(out.rhoresult[13]; digits=2)) kA, rho = $(round(out.rhoresult[1]; digits=3)), flag = $(Int(out.rhoresult[20])), $(out.iend) ray points")
        out.rhoresult[20] == 0 || @warn "beam $(dd.ec_launchers.beam[ibeam].name) exited with flag $(Int(out.rhoresult[20]))"
        push!(beams, Dict(
            "name" => dd.ec_launchers.beam[ibeam].name,
            "intinbeam" => inputs.intinbeam,
            "floatinbeam" => inputs.floatinbeam,
            "ni" => inputs.ni, "nj" => inputs.nj, "npsi" => inputs.npsi,
            "eqdata_sum" => sum(inputs.eqdata), "prdata_sum" => sum(inputs.prdata),
            [string(f) => getfield(out, f) for f in fieldnames(TORBEAM.BeamOutputs)]...,
        ))
    end
    golden = Dict(
        "case" => case,
        "generated" => string(now()),
        "julia" => string(VERSION),
        "libtorbeam" => realpath(TORBEAM.fortran_library()),
        "params" => params_dict,
        "beams" => beams,
    )
    open(joinpath(GOLDENS_DIR, "$case.json"), "w") do io
        JSON.print(io, round10(golden))
    end
end
