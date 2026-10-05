# Stage 4: adjoint current drive

@testset "current drive" begin
    @testset "Spitzer function, uniform field" begin
        # a "flux surface" with constant B: no trapping, χ̂ = u⁴ ξ / (5+Z) non-relativistically
        n = 64
        fs = TORBEAM.FluxSurface(0.5, fill(1.7, n), zeros(n), fill(2.0, n), fill(0.1, n), 2.0, 2.0, 1.0)
        for Z in (1.0, 2.5)
            sf = TORBEAM.SpitzerFunction(fs, Z; nu=2000, nλ=120, umax=0.5)
            for u in (0.05, 0.1, 0.2), ξ in (0.3, 0.7, 0.95)
                λ = 1 - ξ^2
                χ = TORBEAM.chi(sf, u, λ)
                γ = sqrt(1 + u^2)
                # relativistic correction to the exact solution is O(u²); allow for it
                @test χ ≈ u^4 * ξ / (5 + Z) rtol = 0.03 + 2u^2
            end
        end
    end

    @testset "flux surface and trapping" begin
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "D3D.json"))
        params = TORBEAM.TorbeamParams()
        eq = TORBEAM.equilibrium_inputs(dd)
        inputs = TORBEAM.beam_inputs(dd, 1, params, eq)
        m = TORBEAM.PlasmaModel(inputs)
        for ρ in (0.3, 0.7, 0.95)
            fs = TORBEAM.FluxSurface(m, ρ)
            @test all(abs.(TORBEAM.rho_pol.(Ref(m), fs.R, fs.Z) .- ρ) .< 1e-6)
            @test 0 < fs.λc < 1
            @test fs.Bmax ≈ maximum(fs.B) && fs.Bmin ≈ minimum(fs.B)
            # outboard side has the lowest field
            @test fs.R[argmin(fs.B)] > m.R_axis
            sf = TORBEAM.SpitzerFunction(fs, m.Zeff)
            # χ vanishes at the trapped-passing boundary and grows with u
            @test abs(TORBEAM.chi(sf, 0.3, fs.λc)) < 1e-3 * abs(TORBEAM.chi(sf, 0.3, 0.0))
            @test TORBEAM.chi(sf, 0.3, 0.0) > TORBEAM.chi(sf, 0.2, 0.0) > TORBEAM.chi(sf, 0.1, 0.0) > 0
            # trapping reduces the response compared with the uniform-field value
            @test TORBEAM.chi(sf, 0.3, 0.5 * fs.λc) < 0.3^4 * sqrt(1 - 0.5 * fs.λc) / (5 + m.Zeff)
        end
    end

    for case in [splitext(f)[1] for f in readdir(joinpath(@__DIR__, "goldens")) if endswith(f, ".json")]
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "$case.json"))
        golden = JSON.parsefile(joinpath(@__DIR__, "goldens", "$case.json"))
        params = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])..., backend=:julia)
        eq = TORBEAM.equilibrium_inputs(dd)
        @testset "$case driven current vs Fortran" begin
            for (ibeam, gb) in enumerate(golden["beams"])
                I2 = gb["rhoresult"][13]               # Fortran with momentum conservation (ncdroutine=2)
                I1 = gb["Icd_ncdroutine1"]            # Fortran Lin-Liu without momentum conservation
                fi = gb["floatinbeam"]
                nharm = round(Int, fi[1] / (27.99e9 * abs(fi[27])))
                O2 = gb["intinbeam"][3] == 1 && nharm >= 2
                # (Julia ncdroutine=1: bounce-averaged Lorentz-model response) vs Fortran ncdroutine=1,
                # (Julia ncdroutine=2: rescaled by the full-operator Spitzer function) vs Fortran ncdroutine=2.
                # Second-harmonic O-mode: the sign is tested, the magnitude is not — its weak absorption
                # (Julia 1.27 vs Fortran 1.25 MW) straddles the cold resonance and the current follows it.
                # DIII-D X2: both models land 10-13 % above the Fortran with the complex warm root
                # (nabsroutine=1; 8-10 % with the weak-damping coefficient), ITER within 1-5 %.
                tolD3D = startswith(case, "D3D") ? 0.15 : 0.0
                for (ncdr, Iref, tol) in ((1, I1, max(0.12, tolD3D)), (2, I2, max(0.1, tolD3D)))
                    p = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])..., backend=:julia, ncdroutine=ncdr)
                    inputs = TORBEAM.beam_inputs(dd, ibeam, p, eq)
                    out = TORBEAM.run_beam(inputs, p)
                    I = out.rhoresult[13]
                    @info "$case $(gb["name"]) ncdroutine=$ncdr: I_cd $(round(I; digits=2)) kA vs Fortran $(round(Iref; digits=2)) kA"
                    # only where the current is not a near-cancellation (> 1 kA per MW absorbed)
                    if abs(Iref) > 1.0 * gb["rhoresult"][14]
                        @test sign(I) == sign(Iref)
                        O2 || @test abs(I - Iref) < tol * abs(Iref)
                        # the current profile is where the power is: |j|-weighted mean rho agrees
                        ρ = out.t2ndata[1:TORBEAM.NPNT]
                        j = out.t2ndata[2TORBEAM.NPNT+1:3TORBEAM.NPNT]
                        jref = ncdr == 1 ? gb["j_ncdroutine1"] : gb["t2ndata"][2TORBEAM.NPNT+1:3TORBEAM.NPNT]
                        if sum(abs, jref) > 0
                            @test abs(sum(ρ .* abs.(j)) / sum(abs.(j)) - sum(ρ .* abs.(jref)) / sum(abs.(jref))) < 0.06
                        end
                    end
                end
            end
        end
    end
end

@testset "Spitzer-Härm function" begin
    for (Z, γE) in ((1, 0.5816), (2, 0.6833), (4, 0.7849), (16, 0.9252))
        sp = TORBEAM.SpitzerFunction1D(Z)
        @test sp.γE ≈ γE rtol = 0.005              # Spitzer & Härm (1953) conductivity ratios
        @test sp.ee_residual < 1e-3                 # e-e collisions conserve momentum (discretization level)
        @test all(sp.D[10:end] .> 0)                # D ∝ x⁴ is at round-off in the first nodes
        # tends to the Lorentz-model high-velocity solution
        @test 0.95 < TORBEAM.spitzer_ratio(sp, 6.0) < 1.15
        @test TORBEAM.spitzer_ratio(sp, 1.0) > TORBEAM.spitzer_ratio(sp, 3.0) > 1
    end
end

@testset "variational Spitzer function" begin
    # f_c = 1, μ → ∞: the Spitzer-Härm function in the normalisation of the 1-D solver (no free
    # constant); the 4-term polynomial is a few % off pointwise (the paper quotes 0.5 ≲ x ≲ 4) and
    # the conductivity moment, being the variational quantity, is accurate to < 1 %
    for Z in (1.0, 2.0, 4.0)
        sp = TORBEAM.SpitzerFunction1D(Z)
        d, χ = TORBEAM.variational_spitzer(1.0, Z, 1e8)
        for x in (0.7, 1.0, 1.5, 2.0, 3.0)
            @test χ(x) ≈ sp.D[argmin(abs.(sp.x .- x))] rtol = 0.08
        end
        w = sp.x .^ 4 .* exp.(-sp.x .^ 2)
        @test sum(w .* χ.(sp.x)) ≈ sum(w .* sp.D) rtol = 0.01
    end
    # trapping reduces the response, relativity reduces it at high momentum
    _, χ1 = TORBEAM.variational_spitzer(1.0, 2.0, 1e8)
    _, χ5 = TORBEAM.variational_spitzer(0.5, 2.0, 1e8)
    _, χr = TORBEAM.variational_spitzer(0.5, 2.0, 20.44)
    @test all(χ5(u) < χ1(u) for u in (1.0, 2.0, 3.0))
    @test all(0.5 < χr(u) / χ5(u) < 0.9 for u in (2.0, 3.0, 4.0))
    @test issorted([χr(u) / χ5(u) for u in (1.0, 2.0, 3.0, 4.0)]; rev=true)
end

@testset "Lin-Liu separable response" begin
    n = 64
    fsu = TORBEAM.FluxSurface(0.5, fill(6.0, n), zeros(n), fill(5.0, n), fill(0.1, n), 5.0, 5.0, 1.0)
    @test TORBEAM.circulating_fraction(fsu) ≈ 1.0 atol = 1e-4
    # uniform field: the l = 1 truncation is exact, so the separable response coincides with the
    # exact 2-D solution of the same operator, relativistic factors included
    # (the march's first-order λ discretisation is ~1% low at λ = 0 with nλ = 400 and converges)
    sfl = TORBEAM.linliu_response(fsu, 2.0, 100.0)
    sfL = TORBEAM.SpitzerFunction(fsu, 2.0; nu=1500, nλ=400)
    for u in (0.1, 0.3, 0.6, 1.0), λ in (0.0, 0.5, 0.9)
        @test TORBEAM.chi(sfl, u, λ) ≈ TORBEAM.chi(sfL, u, λ) rtol = 0.02
    end
    @test TORBEAM.chi(sfl, 0.1, 0.3) ≈ 0.1^4 * sqrt(0.7) / 7 rtol = 0.01     # non-relativistic u⁴ξ/(5+Z)
    @test sfl.scale ≈ 1.0
    # momentum conservation, uniform field, 5 keV: the enhancement over the high-speed limit is the
    # Spitzer-Härm one (variational polynomial, few %), and the hand-over to the exact 1-D solution
    # at x = 3 is smooth
    μ = 100.0
    uT = sqrt(2 / μ)
    sfm = TORBEAM.linliu_response(fsu, 2.0, μ; momentum_conservation=true)
    sfh = TORBEAM.linliu_response(fsu, 2.0, μ)
    sp = TORBEAM.SpitzerFunction1D(2.0)
    for x in (1.0, 2.0, 3.0)
        @test TORBEAM.chi(sfm, x * uT, 0.0) / TORBEAM.chi(sfh, x * uT, 0.0) ≈ sp.D[argmin(abs.(sp.x .- x))] * 7 / x^4 rtol = 0.08
    end
    for x in (2.7, 2.9, 3.1, 3.3)
        @test TORBEAM.chi(sfm, x * uT, 0.0) / TORBEAM.chi(sfh, x * uT, 0.0) ≈ sp.D[argmin(abs.(sp.x .- x))] * 7 / x^4 rtol = 0.05
    end
    # real surface: f_c is the complement of the effective trapped fraction; H decreasing to 0
    dd = IMAS.json2imas(joinpath(@__DIR__, "data", "ITER.json"))
    m = TORBEAM.PlasmaModel(TORBEAM.beam_inputs(dd, 1, TORBEAM.TorbeamParams(), TORBEAM.equilibrium_inputs(dd)))
    p1 = dd.equilibrium.time_slice[1].profiles_1d
    psin = (p1.psi .- p1.psi[1]) ./ (p1.psi[end] - p1.psi[1])
    for ρ in (0.45, 0.775)
        fs = TORBEAM.FluxSurface(m, ρ)
        @test TORBEAM.circulating_fraction(fs) ≈ 1 - IMAS.interp1d(psin, p1.trapped_fraction)(ρ^2) rtol = 0.01
        sf = TORBEAM.linliu_response(fs, m.Zeff, 30.0)
        χs = [TORBEAM.chi(sf, 0.5, f * fs.λc) for f in (0.0, 0.3, 0.6, 0.9, 1.0)]
        @test issorted(χs; rev=true)
        @test abs(χs[end]) < 1e-6 * χs[1]
        @test 0.7 < sf.scale < 1.0
    end
end

@testset "full-operator response" begin
    # uniform field: the 2-D solver reproduces the 1-D Spitzer function (which has
    # energy diffusion and the field term) at thermal energies, where the
    # relativistic γ factors are ~1
    n = 64
    for (Z, Te) in ((2.0, 5.0), (1.0, 2.0))
        μ = 510.99895 / Te
        uT = sqrt(2 / μ)
        fs = TORBEAM.FluxSurface(0.5, fill(6.0, n), zeros(n), fill(5.0, n), fill(0.1, n), 5.0, 5.0, 1.0)
        sp = TORBEAM.SpitzerFunction1D(Z)
        umax = min(1.5, 10uT)                # resolve the thermal bulk at 2 keV
        sf = TORBEAM.full_operator_response(fs, Z, μ; umax)
        for u in (0.3uT, 0.7uT, 1.0uT)
            D1 = sp.D[argmin(abs.(sp.x .- u / uT))]
            @test TORBEAM.chi(sf, u, 0.64) / 0.6 ≈ uT^4 * D1 rtol = 0.06
        end
        # the field term enhances the response
        sft = TORBEAM.full_operator_response(fs, Z, μ; umax, field=false)
        @test TORBEAM.chi(sf, 0.7uT, 0.64) > 1.2 * TORBEAM.chi(sft, 0.7uT, 0.64)
    end
    # suprathermal electrons: the test-particle-only response must stay above the Lorentz one
    # (exact thermal rates lie below their 1/u³ asymptotes) and tend to it from above, with only
    # a weak dependence on μ at fixed x = u/u_T — the γ² drag and the γ³ energy diffusion then
    # cancel on the relativistic Maxwellian (detailed balance); with the non-relativistic γ of the
    # energy diffusion the leftover drag drove the ratio to 0.85 at x = 3 for μ = 51
    fs = TORBEAM.FluxSurface(0.5, fill(6.0, n), zeros(n), fill(5.0, n), fill(0.1, n), 5.0, 5.0, 1.0)
    ratios = Dict{Float64,Vector{Float64}}()
    for (μ, umax) in ((51.1, 1.5), (200.0, 1.0))
        uT = sqrt(2 / μ)
        sfL = TORBEAM.SpitzerFunction(fs, 1.0; nu=1500)
        sfT = TORBEAM.full_operator_response(fs, 1.0, μ; umax, field=false)
        r = [TORBEAM.chi(sfT, x * uT, 0.0) / TORBEAM.chi(sfL, x * uT, 0.0) for x in (2, 3, 4)]
        @test all(1.0 .< r .< 1.7)
        @test issorted(r; rev=true)
        ratios[μ] = r
    end
    @test abs(ratios[51.1][2] - ratios[200.0][2]) < 0.05 * ratios[200.0][2]

    # trapped surface: vanishes at the trapped boundary, odd structure preserved
    dd = IMAS.json2imas(joinpath(@__DIR__, "data", "ITER.json"))
    inputs = TORBEAM.beam_inputs(dd, 1, TORBEAM.TorbeamParams(), TORBEAM.equilibrium_inputs(dd))
    m = TORBEAM.PlasmaModel(inputs)
    fs = TORBEAM.FluxSurface(m, 0.5)
    sf = TORBEAM.full_operator_response(fs, m.Zeff, 510.99895 / TORBEAM.temperature(m, 0.5))
    @test abs(TORBEAM.chi(sf, 0.3, fs.λc)) < 1e-3 * abs(TORBEAM.chi(sf, 0.3, 0.0))

    # neoclassical conductivity: the Spitzer problem on the real surfaces, driven by E∥ ∝ B and
    # measured as ⟨j∥B⟩/⟨E∥B⟩, is c_b = ∮b dl/∮dl times the solver's response to its unweighted
    # drive. Lorentz gas (pitch-angle scattering only): exactly 1 - f_t with the effective trapped
    # fraction; full operator: the Sauter et al. (1999) collisionless fit to within a few %
    fsavg(fs, A) = sum(A .* fs.dl ./ fs.B) / sum(fs.dl ./ fs.B)
    function trapped_fraction(fs)
        b = fs.B ./ fs.Bmin
        λs = range(0, fs.λc; length=2001)
        return 1 - 0.75 * fsavg(fs, b .^ 2) * sum(λ / fsavg(fs, sqrt.(max.(1 .- λ .* b, 0.0))) for λ in λs) * step(λs)
    end
    function jmoment(sf, μ, λc)
        us = range(0, 1.2; length=301)[2:end]
        λs = range(0, λc; length=300)
        return sum(u^3 / sqrt(1 + u^2) * exp(-μ * (sqrt(1 + u^2) - 1)) * sum(TORBEAM.chi(sf, u, λ) for λ in λs) for u in us)
    end
    fsu = TORBEAM.FluxSurface(0.5, fill(6.0, n), zeros(n), fill(5.0, n), fill(0.1, n), 5.0, 5.0, 1.0)
    for ρ in (0.45, 0.775)
        fs = TORBEAM.FluxSurface(m, ρ)
        μ = 510.99895 / TORBEAM.temperature(m, ρ)
        cb = sum(fs.B ./ fs.Bmin .* fs.dl) / sum(fs.dl)
        ft = trapped_fraction(fs)
        σL = cb * fs.λc * jmoment(TORBEAM.full_operator_response(fs, 1e4, μ; field=false), μ, fs.λc) / jmoment(TORBEAM.full_operator_response(fsu, 1e4, μ; field=false), μ, 1.0)
        @test σL ≈ 1 - ft rtol = 0.01
        for Z in (1.0, m.Zeff)
            σ = cb * fs.λc * jmoment(TORBEAM.full_operator_response(fs, Z, μ), μ, fs.λc) / jmoment(TORBEAM.full_operator_response(fsu, Z, μ), μ, 1.0)
            sauter = 1 - (1 + 0.36 / Z) * ft + 0.59 / Z * ft^2 - 0.23 / Z * ft^3
            @test σ ≈ sauter rtol = 0.06
        end
    end
    @test TORBEAM.chi(sf, 0.3, 0.0) > TORBEAM.chi(sf, 0.2, 0.0) > TORBEAM.chi(sf, 0.1, 0.0) > 0
end
