using TORBEAM
using TORBEAM.IMAS
using JSON
using Test

include("test_model.jl")
include("test_beam_tracing.jl")
include("test_absorption.jl")

const DATA_DIR = joinpath(@__DIR__, "data")
const GOLDENS_DIR = joinpath(@__DIR__, "goldens")

cases = [splitext(f)[1] for f in readdir(GOLDENS_DIR) if endswith(f, ".json")]
@assert !isempty(cases) "no golden data in $GOLDENS_DIR"

if !TORBEAM.fortran_available()
    @warn "TORBEAM Fortran library not found (TORBEAM_DIR): skipping the Fortran regression tests"
end

@testset "TORBEAM" begin
    for case in cases
        dd = IMAS.json2imas(joinpath(DATA_DIR, "$case.json"))
        golden = JSON.parsefile(joinpath(GOLDENS_DIR, "$case.json"))
        params = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])...)
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
            @testset "$case fortran backend" begin
                eq = TORBEAM.equilibrium_inputs(dd)
                for (ibeam, gb) in enumerate(golden["beams"])
                    inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
                    out = TORBEAM.run_beam(inputs, params)
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

            @testset "$case run_torbeam" begin
                TORBEAM.run_torbeam(dd, params)
                @test length(dd.waves.coherent_wave) == nbeam
                for (ibeam, gb) in enumerate(golden["beams"])
                    wv = dd.waves.coherent_wave[ibeam]
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
