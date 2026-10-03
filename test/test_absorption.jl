# Stage 3: absorption, deposition profiles and the Julia backend end to end

import LinearAlgebra: eigvals, Hermitian

@testset "absorption" begin
    @testset "anti-Hermitian tensor" begin
        μ = 510.99895 / 2.0     # 2 keV
        for (X, Y, Nperp, Npar) in ((0.3, 0.52, 0.8, 0.0), (0.3, 0.52, 0.6, 0.3), (0.4, 1.02, 0.7, -0.25), (0.2, 0.3, 0.9, 0.1))
            εa = TORBEAM.antihermitian_tensor(X, Y, Nperp, Npar, μ)
            @test εa ≈ εa' atol = 1e-12 * max(1, maximum(abs, εa))         # Hermitian
            @test all(>=(-1e-12 * max(1, maximum(abs, εa))), eigvals(Hermitian(εa)))  # positive semidefinite
            # isotropic Maxwellian: N∥ -> -N∥ is the reflection z -> -z, under which
            # the polarization flips its z component and the absorption is unchanged
            for e in (TORBEAM.cold_polarization(X, Y, Nperp, Npar), [0, 0, 1 + 0im])
                εa2 = TORBEAM.antihermitian_tensor(X, Y, Nperp, -Npar, μ)
                e2 = [e[1], e[2], -e[3]]
                @test real(e' * εa * e) ≈ real(e2' * εa2 * e2) rtol = 1e-6
            end
        end
        @test iszero(TORBEAM.antihermitian_tensor(0.0, 0.5, 0.8, 0.0, μ))        # no plasma
        @test iszero(TORBEAM.antihermitian_tensor(0.3, 0.3, 0.8, 0.0, μ))        # no resonance (n Y + N∥² < 1 for n ≤ 3)
        # cold polarization is a null vector of D and |e| = 1
        for (X, Y, Nperp, Npar, mode) in ((0.3, 0.52, 0.8, 0.0, 1), (0.3, 0.52, 0.8, 0.0, -1))
            n2 = TORBEAM.refractive_index2(X, Y, Npar^2 / (Npar^2 + Nperp^2), mode)
            # rescale Nperp so (Nperp, Npar) is on the dispersion surface
            Np = sqrt(max(n2 - Npar^2, 0))
            e = TORBEAM.cold_polarization(X, Y, Np, Npar)
            S, D, P = TORBEAM.cold_tensor(X, Y)
            N2 = Np^2 + Npar^2
            M = [S-Npar^2 -im*D Np*Npar; im*D S-N2 0; Np*Npar 0 P-Np^2]
            @test norm(M * e) < 1e-10
            @test norm(e) ≈ 1
            mode == 1 && @test abs(e[3]) > 0.9      # O-mode: E along B at perpendicular propagation
            mode == -1 && @test abs(e[3]) < 1e-10   # X-mode: E ⟂ B
        end
    end

    for case in [splitext(f)[1] for f in readdir(joinpath(@__DIR__, "goldens")) if endswith(f, ".json")]
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "$case.json"))
        golden = JSON.parsefile(joinpath(@__DIR__, "goldens", "$case.json"))
        params = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])..., backend=:julia)
        nprofv = params.nprofv
        eq = TORBEAM.equilibrium_inputs(dd)
        @testset "$case julia backend vs Fortran" begin
            for (ibeam, gb) in enumerate(golden["beams"])
                inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
                out = TORBEAM.run_beam(inputs, params)
                Pin = inputs.floatinbeam[24]
                @test out.rhoresult[20] == 0
                # total absorbed power
                @test out.rhoresult[14] ≈ gb["rhoresult"][14] atol = 0.12 * Pin
                # R, Z [cm] of the absorption maximum along the central ray
                @test hypot(out.rhoresult[2] - gb["rhoresult"][2], out.rhoresult[3] - gb["rhoresult"][3]) < 12
                # deposition profile: integrates to the absorbed power on our volumes,
                # and its median in rho agrees with the Fortran's (the dP/dV maximum is
                # not a robust location near the axis, where dV/dρ -> 0)
                function cumulative(t2n, vp)
                    ρ = t2n[1:TORBEAM.NPNT]
                    dPdV = t2n[TORBEAM.NPNT+1:2TORBEAM.NPNT]
                    V = TORBEAM.cubic_resample(vp[1:nprofv], vp[nprofv+1:2nprofv], ρ)
                    dV = [k == 1 ? V[2] - V[1] : V[k] - V[k-1] for k in eachindex(V)]
                    return ρ, dPdV, cumsum(dPdV .* dV)
                end
                ρ, dPdV, c = cumulative(out.t2ndata, out.volprof)
                _, _, cg = cumulative(gb["t2ndata"], gb["volprof"])
                @test c[end] ≈ out.rhoresult[14] rtol = 0.05
                @test all(>=(0), dPdV)
                median(c) = ρ[findfirst(>=(0.5 * c[end]), c)]
                @test abs(median(c) - median(cg)) < 0.05
                @info "$case $(gb["name"]): P_abs $(round(out.rhoresult[14]; digits=3)) vs $(round(gb["rhoresult"][14]; digits=3)) MW, median rho $(round(median(c); digits=3)) vs $(round(median(cg); digits=3)), (R,Z) of max ($(round(out.rhoresult[2]; digits=1)), $(round(out.rhoresult[3]; digits=1))) vs ($(round(gb["rhoresult"][2]; digits=1)), $(round(gb["rhoresult"][3]; digits=1))) cm"
                if ibeam == 1
                    # flux-surface volumes against the Fortran's (its rho grid extends past 1)
                    for ρk in (0.3, 0.6, 0.9, 1.0)
                        Vj = TORBEAM.cubic_resample(out.volprof[1:nprofv], out.volprof[nprofv+1:2nprofv], [ρk])[1]
                        Vf = TORBEAM.cubic_resample(gb["volprof"][1:nprofv], gb["volprof"][nprofv+1:2nprofv], [ρk])[1]
                        @test Vj ≈ Vf rtol = 0.03
                    end
                end
            end
        end
        @testset "$case run_torbeam with the julia backend" begin
            n = TORBEAM.run_torbeam(dd, params)
            @test n == length(dd.ec_launchers.beam)
            wv = dd.waves.coherent_wave[1]
            @test wv.global_quantities[].power ≈ 1e6 * golden["beams"][1]["rhoresult"][14] rtol = 0.12
            @test length(wv.beam_tracing[].beam) == params.n_ray
            @test length(findall(:ec, dd.core_sources.source)) == n
        end
    end
end
