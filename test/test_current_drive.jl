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
                # Second-harmonic O-mode is skipped: its weak absorption straddles the cold resonance,
                # and which side dominates (hence the sign) depends on the warm corrections of stage 3b.
                for (ncdr, Iref, tol) in ((1, I1, 0.35), (2, I2, 0.65))
                    p = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])..., backend=:julia, ncdroutine=ncdr)
                    inputs = TORBEAM.beam_inputs(dd, ibeam, p, eq)
                    out = TORBEAM.run_beam(inputs, p)
                    I = out.rhoresult[13]
                    @info "$case $(gb["name"]) ncdroutine=$ncdr: I_cd $(round(I; digits=2)) kA vs Fortran $(round(Iref; digits=2)) kA"
                    # only where the current is not a near-cancellation (> 1 kA per MW absorbed)
                    if abs(Iref) > 1.0 * gb["rhoresult"][14] && !O2
                        @test sign(I) == sign(Iref)
                        @test abs(I - Iref) < tol * abs(Iref)
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
