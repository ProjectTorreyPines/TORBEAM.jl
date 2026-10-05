# Compare the Julia backend with the Fortran library on the golden beams.
#
# For every beam of test/goldens/<case>.json this prints the Fortran result,
# the Julia result with the switches that select the same reduced models
# (`ncdroutine` 1 and 2, `nabsroutine` 1, `nprofcalc` 1, the defaults) and the
# Julia result with the higher-fidelity options (`ncdroutine` 3 and 4,
# `nprofcalc` 2). The first comparison is the one to check when deciding
# whether to trust the backend; the second shows what the better-founded
# solvers change and why (see the switch table in the README).
#
# Needs only the JSON files (no Fortran library, no FUSE). It runs every beam
# four times, 20-40 min on one node:
#
#   julia --project=<env with JSON and this package dev'ed> test/goldens/compare.jl [D3D] [ITER]
#
# Set COMPARE_HIFI=0 to skip the higher-fidelity runs.

using TORBEAM
using TORBEAM.IMAS
using JSON
using Printf

const NPNT = TORBEAM.NPNT
const HIFI = get(ENV, "COMPARE_HIFI", "1") != "0"

"""
    profile_stats(ρ, dPdV, V)

Median ρ, 16-84 % width and peak of a deposition profile on the flux-volume grid `V(ρ)`
"""
function profile_stats(ρ, dPdV, V)
    dV = [k == 1 ? V[1] : V[k] - V[k-1] for k in 1:NPNT]
    c = cumsum(dPdV .* dV)
    tot = c[end]
    tot > 0 || return (NaN, NaN, 0.0)
    q(f) = ρ[findfirst(>=(f * tot), c)]
    return q(0.5), q(0.84) - q(0.16), maximum(dPdV)
end

function stats(out, nprofv)
    ρ = out.t2ndata[1:NPNT]
    V = TORBEAM.cubic_resample(out.volprof[1:nprofv], out.volprof[nprofv+1:2nprofv], ρ)
    med, w, pk = profile_stats(ρ, out.t2ndata[NPNT+1:2NPNT], V)
    return (P=out.rhoresult[14], med, w, pk, I=out.rhoresult[13])
end

# the goldens are plain Dicts with the same layout
function stats(gb::Dict, nprofv)
    ρ = Float64.(gb["t2ndata"][1:NPNT])
    V = TORBEAM.cubic_resample(Float64.(gb["volprof"][1:nprofv]), Float64.(gb["volprof"][nprofv+1:2nprofv]), ρ)
    med, w, pk = profile_stats(ρ, Float64.(gb["t2ndata"][NPNT+1:2NPNT]), V)
    return (P=gb["rhoresult"][14], med, w, pk, I=gb["rhoresult"][13])
end

fmt(x, d) = isnan(x) ? "-" : string(round(x; digits=d))
pair(a, b, d) = fmt(a, d) * " / " * fmt(b, d)

"""
    run_variant(dd, eq, golden, ibeams; kw...)

Run every beam in `ibeams` with `TorbeamParams(; backend=:julia, kw...)` on top of the golden
parameters, sharing one `BackendCache`; returns the stats and the wall time per beam
"""
function run_variant(dd, eq, golden, ibeams; kw...)
    p = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])..., backend=:julia, kw...)
    inputs = [TORBEAM.beam_inputs(dd, i, p, eq) for i in ibeams]
    cache = TORBEAM.backend_cache(inputs[1], p)
    return map(inputs) do inp
        t = @elapsed out = TORBEAM.run_beam(inp, p; cache)
        merge(stats(out, p.nprofv), (t=t,))
    end
end

cases = isempty(ARGS) ? ["D3D", "ITER"] : ARGS
for case in cases
    dd = IMAS.json2imas(joinpath(@__DIR__, "..", "data", "$case.json"))
    golden = JSON.parsefile(joinpath(@__DIR__, "$case.json"))
    eq = TORBEAM.equilibrium_inputs(dd)
    nprofv = golden["params"]["nprofv"]
    beams = golden["beams"]
    ibeams = eachindex(beams)
    ref = [stats(gb, nprofv) for gb in beams]
    I1ref = [gb["Icd_ncdroutine1"] for gb in beams]

    println("\n## $case: Fortran vs Julia with the same reduced models (ncdroutine 1/2, nabsroutine 1, nprofcalc 1)\n")
    j1 = run_variant(dd, eq, golden, ibeams; ncdroutine=1)
    j2 = run_variant(dd, eq, golden, ibeams; ncdroutine=2)
    println("| beam | P_abs [MW] | median ρ | width 16-84 % | peak dP/dV [MW/m³] | I ncdroutine=1 [kA] | I ncdroutine=2 [kA] | s/beam |")
    println("|---|---|---|---|---|---|---|---|")
    for (k, gb) in enumerate(beams)
        r, a, b = ref[k], j1[k], j2[k]
        println("| $(gb["name"]) | $(pair(r.P, b.P, 3)) | $(pair(r.med, b.med, 3)) | $(pair(r.w, b.w, 3)) | $(pair(r.pk, b.pk, 3)) | $(pair(I1ref[k], a.I, 2)) | $(pair(r.I, b.I, 2)) | $(round(b.t; digits=1)) |")
        flush(stdout)
    end
    println("\n(Fortran / Julia; the Fortran's ncdroutine=1 current is stored separately in the golden.)")

    HIFI || continue
    println("\n## $case: Julia higher-fidelity options (nprofcalc 2, ncdroutine 3 exact 2-D Lorentz / 4 full linearized operator)\n")
    j3 = run_variant(dd, eq, golden, ibeams; ncdroutine=3, nprofcalc=2)
    j4 = run_variant(dd, eq, golden, ibeams; ncdroutine=4, nprofcalc=2)
    println("| beam | median ρ | width 16-84 % | peak dP/dV [MW/m³] | I ncdroutine=3 vs Fortran ncdroutine=1 [kA] | I ncdroutine=4 vs Fortran ncdroutine=2 [kA] |")
    println("|---|---|---|---|---|---|")
    for (k, gb) in enumerate(beams)
        r, a, b = ref[k], j3[k], j4[k]
        println("| $(gb["name"]) | $(pair(r.med, b.med, 3)) | $(pair(r.w, b.w, 3)) | $(pair(r.pk, b.pk, 3)) | $(pair(I1ref[k], a.I, 2)) | $(pair(r.I, b.I, 2)) |")
        flush(stdout)
    end
    println("""
    (Fortran / Julia. nprofcalc=2 spreads each absorption step on the local resonance surface
    instead of the vertical plane. ncdroutine=3 solves the bounce-averaged adjoint equation
    exactly in (u, λ) with the high-velocity operator of the Lin-Liu model and no momentum
    conservation, so it is compared with the Fortran's ncdroutine=1 current. ncdroutine=4 adds
    the thermal collision rates, energy diffusion and the e-e field term (validated against the
    Spitzer-Härm and Sauter conductivities) and is compared with the Fortran's momentum-
    conserving ncdroutine=2. The README's switch table lists the typical differences.)
    """)
end
