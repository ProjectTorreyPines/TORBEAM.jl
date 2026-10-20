# Regression tests against golden data produced by the TORBEAM Fortran library
# (test/goldens/<case>.json) for the trimmed `dd`s in test/data/<case>.json.
#
# - The Julia-side assembly of the TORBEAM inputs is always checked.
# - The Fortran backend and the end-to-end `run_torbeam` are checked when the
#   library can be found (`TORBEAM_DIR`; on omega: `module load torbeam/gcc11.x`),
#   and the wall time per beam is printed next to the time recorded in the golden.
#   These tests are meant to be run on omega and are deliberately not part of CI.

using TORBEAM
using TORBEAM.IMAS
using JSON
using Test
using Printf

const DATA_DIR = joinpath(@__DIR__, "data")
const GOLDENS_DIR = joinpath(@__DIR__, "goldens")

cases = sort([splitext(f)[1] for f in readdir(GOLDENS_DIR) if endswith(f, ".json")])
@assert !isempty(cases) "no golden data in $GOLDENS_DIR"

if !TORBEAM.fortran_available()
    @warn "TORBEAM Fortran library not found (TORBEAM_DIR): skipping the Fortran regression tests"
end

"""
    golden_params(golden)

`TorbeamParams` with the values stored in the golden (unknown keys are ignored)
"""
function golden_params(golden::AbstractDict)
    kw = Dict{Symbol,Any}()
    for (k, v) in golden["params"]
        f = Symbol(k)
        f in fieldnames(TORBEAM.TorbeamParams) || continue
        kw[f] = v isa String ? Symbol(v) : v
    end
    return TORBEAM.TorbeamParams(; kw...)
end

"""
    print_timing(case, golden, timings)

Table of the measured wall time per beam next to the time recorded in the golden
"""
function print_timing(case::String, golden::AbstractDict, timings::Vector)
    println("\nFortran wall time per beam, $case (golden recorded on $(get(golden, "host", "?")), $(golden["generated"]))")
    @printf("%-28s %10s %10s %7s\n", "beam", "golden [s]", "now [s]", "ratio")
    tg = tn = 0.0
    for (name, t_golden, t_now) in timings
        tg += t_golden
        tn += t_now
        @printf("%-28s %10.2f %10.2f %7.2f\n", name, t_golden, t_now, t_now / t_golden)
    end
    @printf("%-28s %10.2f %10.2f %7.2f\n", "total", tg, tn, tn / tg)
    return nothing
end

@testset "TORBEAM" begin
    for case in cases
        dd = IMAS.json2imas(joinpath(DATA_DIR, "$case.json"))
        golden = JSON.parsefile(joinpath(GOLDENS_DIR, "$case.json"))
        params = golden_params(golden)
        nbeam = length(dd.ec_launchers.beam)
        @test nbeam == length(golden["beams"])

        @testset "$case inputs" begin
            # the Julia-side assembly of the TORBEAM inputs must not drift
            eq = TORBEAM.equilibrium_inputs(dd)
            for (ibeam, gb) in enumerate(golden["beams"])
                inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
                @test inputs.ni == gb["ni"] && inputs.nj == gb["nj"] && inputs.npsi == gb["npsi"]
                @test inputs.intinbeam == gb["intinbeam"]
                @test inputs.floatinbeam ≈ gb["floatinbeam"] rtol = 1e-9  # goldens are rounded to 10 significant digits
                @test sum(inputs.eqdata) ≈ gb["eqdata_sum"] rtol = 1e-9
                @test sum(inputs.prdata) ≈ gb["prdata_sum"] rtol = 1e-9
            end
        end

        if TORBEAM.fortran_available()
            timings = []
            @testset "$case fortran backend" begin
                eq = TORBEAM.equilibrium_inputs(dd)
                for (ibeam, gb) in enumerate(golden["beams"])
                    inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
                    t = @elapsed out = TORBEAM.run_beam(inputs, params)
                    push!(timings, (gb["name"], Float64(get(gb, "time_s", NaN)), t))
                    @test out.iend == gb["iend"]
                    @test out.rhoresult ≈ gb["rhoresult"] rtol = 1e-6
                    @test out.t1data ≈ gb["t1data"] rtol = 1e-6
                    @test out.t1tdata ≈ gb["t1tdata"] rtol = 1e-6
                    @test out.t2ndata ≈ gb["t2ndata"] rtol = 1e-6
                    @test out.volprof ≈ gb["volprof"] rtol = 1e-6
                    traj = TORBEAM.ray_trajectories(out)
                    @test size(traj) == (15, min(out.iend, TORBEAM.NTRAJ))
                end
            end
            print_timing(case, golden, timings)

            @testset "$case run_torbeam" begin
                TORBEAM.run_torbeam(dd, params)
                @test length(dd.waves.coherent_wave) == nbeam
                for (ibeam, gb) in enumerate(golden["beams"])
                    wv = dd.waves.coherent_wave[ibeam]
                    # rhoresult(13) is the absorbed power [MW] and rhoresult(12) the driven current [kA] (0-based)
                    @test wv.global_quantities[].power ≈ 1e6 * gb["rhoresult"][14] rtol = 1e-6
                    @test wv.global_quantities[].current_tor ≈ 1e3 * gb["rhoresult"][13] rtol = 1e-6
                    @test length(wv.beam_tracing[].beam) == params.n_ray
                    @test length(wv.beam_tracing[].beam[1].position.r) == min(gb["iend"], TORBEAM.NTRAJ)
                end
                @test length(findall(:ec, dd.core_sources.source)) == nbeam
            end
        end
    end
end
