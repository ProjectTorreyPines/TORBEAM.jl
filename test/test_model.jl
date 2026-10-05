# The plasma/equilibrium model layer reproduces its inputs and is differentiable

import ForwardDiff
import Interpolations

@testset "PlasmaModel" begin
    for case in [splitext(f)[1] for f in readdir(joinpath(@__DIR__, "data")) if endswith(f, ".json")]
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "$case.json"))
        params = TORBEAM.TorbeamParams()
        eq = TORBEAM.equilibrium_inputs(dd)
        inputs = TORBEAM.beam_inputs(dd, 1, params, eq)
        m = TORBEAM.PlasmaModel(inputs)

        eqt = dd.equilibrium.time_slice[]
        eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
        gq = eqt.global_quantities
        cp1d = dd.core_profiles.profiles_1d[]

        @testset "$case grid reproduced" begin
            # cubic splines interpolate the nodes exactly
            for i in 1:7:length(m.R), j in 1:5:length(m.Z)
                R, Z = m.R[i], m.Z[j]
                @test TORBEAM.psi(m, R, Z) ≈ eqt2d.psi[i, j] rtol = 1e-10
                BR, Bφ, BZ = TORBEAM.B_cyl(m, R, Z)
                @test BR ≈ eqt2d.b_field_r[i, j] atol = 1e-10
                @test Bφ ≈ eqt2d.b_field_tor[i, j] rtol = 1e-10
                @test BZ ≈ eqt2d.b_field_z[i, j] atol = 1e-10
            end
            # between nodes: close to what IMAS's own psi spline gives
            _, _, PSI = IMAS.ψ_interpolant(eqt2d)
            for (R, Z) in ((0.5 * (m.R[10] + m.R[11]), 0.5 * (m.Z[20] + m.Z[21])), (m.R_axis + 0.3, m.Z_axis - 0.2))
                @test TORBEAM.psi(m, R, Z) ≈ PSI(R, Z) rtol = 1e-6
            end
        end

        @testset "$case magnetic axis" begin
            @test abs(m.R_axis - gq.magnetic_axis.r) < 0.01
            @test abs(m.Z_axis - gq.magnetic_axis.z) < 0.01
            @test abs(m.ψ_axis - gq.psi_axis) < 1e-3 * abs(gq.psi_boundary - gq.psi_axis)
            @test m.ψ_boundary == gq.psi_boundary
            @test TORBEAM.rho_pol(m, m.R_axis, m.Z_axis) == 0
            @test TORBEAM.rho_pol(m, m.R_axis + 0.9 * m.a, m.Z_axis) > 0.5
        end

        @testset "$case profiles" begin
            ψ = cp1d.grid.psi
            ρ = sqrt.((ψ .- ψ[1]) ./ (ψ[end] .- ψ[1]))
            for k in 1:10:length(ρ)
                @test TORBEAM.density(m, ρ[k]) ≈ cp1d.electrons.density[k] rtol = 1e-3
                @test TORBEAM.temperature(m, ρ[k]) ≈ cp1d.electrons.temperature[k] * 1e-3 rtol = 1e-3
            end
            @test m.ρ_edge ≈ 1.0 atol = 1e-12
            # smooth decay into the vacuum region
            @test TORBEAM.density(m, 1.0) ≈ cp1d.electrons.density[end] rtol = 1e-6
            @test TORBEAM.density(m, 1.0) > TORBEAM.density(m, 1.05) > TORBEAM.density(m, 1.2) > 0
            @test TORBEAM.density(m, 1.5) < 1e-6 * TORBEAM.density(m, 1.0)
            @test TORBEAM.temperature(m, 1.3) ≈ TORBEAM.temperature(m, 1.0) rtol = 1e-12
            @test m.Zeff ≈ cp1d.zeff[1]
            @test m.B0 ≈ gq.vacuum_toroidal_field.b0
        end

        @testset "$case cartesian and derivatives" begin
            R, Z = m.R_axis + 0.4 * m.a, m.Z_axis + 0.1
            for φ in (0.0, 0.7, -2.0)
                x, y = R * cos(φ), R * sin(φ)
                s = TORBEAM.state(m, x, y, Z)
                BR, Bφ, BZ = TORBEAM.B_cyl(m, R, Z)
                @test s.Bmag ≈ sqrt(BR^2 + Bφ^2 + BZ^2) rtol = 1e-12
                @test s.R ≈ R && s.Z == Z
                @test s.ψn ≈ TORBEAM.rho_pol(m, R, Z)^2
                @test s.ne ≈ TORBEAM.density_ψn(m, s.ψn) && s.Te ≈ TORBEAM.temperature_ψn(m, s.ψn)
                # the field is invariant under toroidal rotation
                @test hypot(s.B[1], s.B[2]) ≈ hypot(BR, Bφ) rtol = 1e-12
                # ForwardDiff goes through the whole state evaluation
                g = ForwardDiff.gradient(p -> TORBEAM.state(m, p[1], p[2], p[3]).ne, [x, y, Z])
                @test all(isfinite, g)
                gψ = ForwardDiff.gradient(p -> TORBEAM.psi(m, p[1], p[2]), [R, Z])
                @test gψ ≈ Interpolations.gradient(m.ψ, R, Z) rtol = 1e-10
            end
        end
    end
end
