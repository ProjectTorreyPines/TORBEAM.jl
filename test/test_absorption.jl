# Absorption, deposition profiles and the Julia backend end to end

import LinearAlgebra: eigvals, Hermitian, dot, norm, I

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

    @testset "harmonic vector" begin
        # the harmonic matrix is the outer product of the gyro-phase Fourier coefficient of the
        # velocity, w = ∫ (cos φ, sin φ, u∥/u⊥) e^{i(nφ - b sin φ)} dφ / 2π, evaluated here by
        # quadrature: this pins the relative signs of all three components (the cold limit only
        # sees the x–y block; the x–z and y–z entries are odd in u∥ and set the ECCD direction)
        φ = range(0, 2π; length=513)[1:end-1]
        for (n, b, upar, uperp) in ((1, 0.3, 0.4, 0.5), (2, 1.2, -0.3, 0.7), (3, 0.05, 0.6, 0.2), (-1, 0.8, 0.2, 0.4))
            wnum = [sum(f.(φ) .* exp.(im .* (n .* φ .- b .* sin.(φ)))) / length(φ) for f in (cos, sin, _ -> upar / uperp)]
            @test maximum(abs, wnum * wnum' - TORBEAM.bessel_matrix(n, b, upar, uperp)) < 1e-10
            @test maximum(abs, wnum - TORBEAM.harmonic_vector(n, b, upar, uperp)) < 1e-10
        end
        @test TORBEAM.harmonic_vector(1, 0.0, 0.3, 0.5) ≈ [0.5, 0.5im, 0.0]
    end

    @testset "complex warm root" begin
        # the analytic continuation of the tensors in N⊥ reduces to the real-argument values
        μ = 510.99895 / 2.0
        for f in (TORBEAM.antihermitian_tensor, TORBEAM.hermitian_tensor)
            @test f(0.35, 0.505, complex(0.64, 0.0), 0.1, μ) ≈ f(0.35, 0.505, 0.64, 0.1, μ) atol = 1e-12
        end
        # weakly damped point (fundamental O-mode at 2 keV, far from the cold resonance): the
        # complex root's imaginary part agrees with the weak-damping value to first order
        X, Y, Npar, mode = 0.3, 0.97, 0.3, 1
        n2 = TORBEAM.refractive_index2(X, Y, Npar^2 / (Npar^2 + 0.5), mode)
        Np = sqrt(max(n2 - Npar^2, 1e-6))
        for _ in 1:30
            n2 = TORBEAM.refractive_index2(X, Y, Npar^2 / (Npar^2 + Np^2), mode)
            Np = sqrt(max(n2 - Npar^2, 1e-6))
        end
        Npw, Dw = TORBEAM.warm_dispersion(X, Y, Npar, μ, mode, Np)
        e = TORBEAM.null_vector(Dw)
        num = real(dot(e, TORBEAM.antihermitian_tensor(X, Y, Npw, Npar, μ) * e))
        λ(Nq) = (Nv = [Nq, 0.0, Npar]; real(dot(e, (Nv * Nv' - (Nq^2 + Npar^2) * I + TORBEAM.hermitian_tensor(X, Y, Nq, Npar, μ)) * e)))
        h = 1e-4
        κ0 = -num / ((λ(Npw + h) - λ(Npw - h)) / (2h))
        @test 0 < κ0 < 0.01 * Npw
        z = TORBEAM.complex_warm_root(X, Y, Npar, μ, mode, Npw, κ0)
        @test z !== nothing
        @test imag(z) ≈ κ0 rtol = 0.1
        @test real(z) ≈ Npw rtol = 0.01
        # strongly damped: the root moves off the weak-damping value but stays a root
        z2 = TORBEAM.complex_warm_root(0.35, 0.505, 0.1, μ, -1, 0.6469, 0.2556)
        @test z2 !== nothing && imag(z2) > 0.1
    end

    @testset "warm dispersion" begin
        # Hermitian tensor: cold limit, and against a brute-force complex-shift evaluation
        for (X, Y, Nperp, Npar) in ((0.3, 0.52, 0.8, 0.0), (0.3, 0.49, 0.6, 0.3), (0.2, 0.3, 0.9, 0.1))
            εh = TORBEAM.hermitian_tensor(X, Y, Nperp, Npar, 510.99895 / 0.01)
            S, D, P = TORBEAM.cold_tensor(X, Y)
            @test maximum(abs, εh - [S -im*D 0; im*D S 0; 0 0 P]) < 1e-3
            @test εh ≈ εh' atol = 1e-12
        end
        function brute(X, Y, Nperp, Npar, μ, δ)
            ε = Matrix{ComplexF64}(TORBEAM.I, 3, 3)
            fnorm = μ / (4π * TORBEAM.besselkx(2, μ))
            L = min(sqrt((1 + 32 / μ)^2 - 1), 3.0)
            for n in -1:3
                Iperp, _ = TORBEAM.quadgk(0.0, L; rtol=1e-6) do uperp
                    b = Nperp * uperp / Y
                    Jn = TORBEAM.besselj(n, b)
                    Jnp = 0.5 * (TORBEAM.besselj(n - 1, b) - TORBEAM.besselj(n + 1, b))
                    nJb = b > 1e-12 ? n * Jn / b : (abs(n) == 1 ? 0.5 * sign(n) : 0.0)
                    Iu, _ = TORBEAM.quadgk(-L, L; rtol=1e-7, maxevals=200000) do upar
                        γ = sqrt(1 + uperp^2 + upar^2)
                        w = [nJb, im * Jnp, Jn * upar / uperp]
                        (uperp^2 * (-μ * fnorm * exp(-μ * (γ - 1))) / γ) .* (w * w') ./ (γ - n * Y - Npar * upar + im * δ)
                    end
                    2π * uperp .* Iu
                end
                ε .+= X .* Iperp
            end
            return ε
        end
        # the complex shift biases the integral linearly in δ: extrapolate 2ε(δ) - ε(2δ);
        # the remaining Hermitian difference at 20 keV is the brute's u cutoff
        for (X, Y, Nperp, Npar, Te, tolh) in ((0.35, 0.505, 0.64, 0.1, 2.0, 0.03), (0.3, 1.01, 0.83, 0.2, 20.0, 0.06))
            μ = 510.99895 / Te
            εh = TORBEAM.hermitian_tensor(X, Y, Nperp, Npar, μ)
            εa = TORBEAM.antihermitian_tensor(X, Y, Nperp, Npar, μ)
            εb = 2 * brute(X, Y, Nperp, Npar, μ, 1e-3) - brute(X, Y, Nperp, Npar, μ, 2e-3)
            @test maximum(abs, (εb + εb') / 2 - εh) < tolh * maximum(abs, εh - TORBEAM.I)
            @test maximum(abs, (εb - εb') / 2im - εa) < 0.02 * maximum(abs, εa)
        end
        # warm N⊥: close to the cold root far from the resonance, shifts by several % near it
        μ = 510.99895 / 2.0
        for (Y, shift) in ((0.3, false), (0.49, true), (0.52, true))
            X, Npar, mode = 0.35, 0.1, -1
            n2 = 1.0
            Np = 0.6
            for _ in 1:40
                n2 = TORBEAM.refractive_index2(X, Y, Npar^2 / (Npar^2 + Np^2), mode)
                Np = sqrt(max(n2 - Npar^2, 1e-6))
            end
            Npw, Dw = TORBEAM.warm_dispersion(X, Y, Npar, μ, mode, Np)
            e = TORBEAM.null_vector(Dw)
            @test norm(Dw * e) < 1e-6 * norm(Dw)
            @test shift ? abs(Npw - Np) > 0.02 : abs(Npw - Np) < 5e-3
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
                # smooth profile: rms relative second difference (10-bin stride) over the bins
                # above 30 % of the peak
                kpk = argmax(dPdV)
                idx = [i for i in 11:TORBEAM.NPNT-10 if dPdV[i] > 0.3 * dPdV[kpk]]
                rough = sqrt(sum(((dPdV[i-10] - 2dPdV[i] + dPdV[i+10]) / dPdV[i])^2 for i in idx) / length(idx))
                @test rough < 0.01
                # the maximum of dP/dV lies where the Fortran's is (no spike at the axis)
                gdPdV = gb["t2ndata"][TORBEAM.NPNT+1:2TORBEAM.NPNT]
                @test abs(ρ[kpk] - ρ[argmax(gdPdV)]) < 0.04
                quantile(c, f) = ρ[findfirst(>=(f * c[end]), c)]
                median(c) = quantile(c, 0.5)
                @test abs(median(c) - median(cg)) < 0.02
                # 16-84 % width and peak of dP/dV against the Fortran's
                @test quantile(c, 0.84) - quantile(c, 0.16) ≈ quantile(cg, 0.84) - quantile(cg, 0.16) rtol = 0.15
                @test maximum(dPdV) ≈ maximum(gdPdV) rtol = 0.2
                # total absorbed power
                @test out.rhoresult[14] ≈ gb["rhoresult"][14] rtol = 0.03
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
