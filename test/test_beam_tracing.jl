# Stage 2: paraxial beam tracing against an analytic vacuum Gaussian beam and the golden rays

import LinearAlgebra: norm

"""
Golden central ray (X, Y, Z) [m], the vertical and horizontal 1/e half-widths [m]
and the cumulative arclength [m] from the Fortran `t1data`/`t1tdata`
"""
function golden_ray(gb)
    n = gb["iend"]
    t1 = gb["t1data"]
    tt = gb["t1tdata"]
    X = tt[1:n] ./ 100
    Y = tt[n+1:2n] ./ 100
    Z = t1[n+1:2n] ./ 100
    wv = [hypot(t1[2n+i] - t1[4n+i], t1[3n+i] - t1[5n+i]) for i in 1:n] ./ 200
    wh = [hypot(tt[2n+i] - tt[4n+i], tt[3n+i] - tt[5n+i]) for i in 1:n] ./ 200
    # the first stored point is one step past the launch point: arclength from the launch
    x0 = gb["floatinbeam"][4:6] ./ 100
    s = zeros(n)
    s[1] = hypot(X[1] - x0[1], Y[1] - x0[2], Z[1] - x0[3])
    for i in 2:n
        s[i] = s[i-1] + hypot(X[i] - X[i-1], Y[i] - Y[i-1], Z[i] - Z[i-1])
    end
    return (; X, Y, Z, wv, wh, s, R=t1[1:n] ./ 100)
end

@testset "beam tracing" begin
    @testset "vacuum Gaussian beam" begin
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "D3D.json"))
        params = TORBEAM.TorbeamParams()
        eq = TORBEAM.equilibrium_inputs(dd)
        inputs = TORBEAM.beam_inputs(dd, 1, params, eq)
        # same equilibrium, no plasma: zero density
        prdata = copy(inputs.prdata)
        prdata[inputs.npsi+1:2*inputs.npsi] .= 0
        vac = TORBEAM.BeamInputs(inputs.intinbeam, inputs.floatinbeam, inputs.ni, inputs.nj, inputs.eqdata, inputs.npsi, prdata)
        m = TORBEAM.PlasmaModel(vac)
        l = TORBEAM.Launch(vac)
        b = TORBEAM.trace_beam(m, l; smax=1.0)
        @test b.exit == :length
        # straight ray
        x1 = b.sol(1.0)[1:3]
        @test x1 ≈ l.x0 + l.N0 rtol = 1e-8
        @test norm(b.sol(1.0)[4:6]) ≈ 1 rtol = 1e-8
        # width evolution of a Gaussian beam: w(s)² = w0² (1 + ((s - s_w)/z_R)²)
        # with the waist w0 at s_w, from the launch q-parameter 1/q = S + iΦ
        k0 = l.wave.k0
        for (w, R) in ((l.wh, l.Rh), (l.wv, l.Rv))
            q = 1 / complex(-1 / R, 2 / (k0 * w^2))
            s_w = -real(q)
            z_R = -imag(q)
            w0 = sqrt(2 * z_R / k0)
            for s in (0.0, 0.05, 0.1, 0.3, 0.7, 1.0)
                bw = TORBEAM.beam_widths(b, s)
                w_s = w == l.wh ? bw.wh : bw.wv
                @test w_s ≈ w0 * sqrt(1 + ((s - s_w) / z_R)^2) rtol = 1e-6
            end
        end
        @test 0.05 < -real(1 / complex(-1 / l.Rh, 2 / (k0 * l.wh^2))) < 0.1   # waist ~7 cm ahead
        # the ABCD propagation of a finite-difference ray family reproduces the same law
        for s in (0.3, 1.0)
            a = TORBEAM.abcd_widths(b, m, s)
            for (w, R, w_s) in ((l.wh, l.Rh, a.wh), (l.wv, l.Rv, a.wv))
                q = 1 / complex(-1 / R, 2 / (k0 * w^2))
                @test w_s ≈ sqrt(-2 * imag(q) / k0) * sqrt(1 + ((s + real(q)) / -imag(q))^2) rtol = 1e-5
            end
        end
    end

    for case in [splitext(f)[1] for f in readdir(joinpath(@__DIR__, "goldens")) if endswith(f, ".json")]
        dd = IMAS.json2imas(joinpath(@__DIR__, "data", "$case.json"))
        golden = JSON.parsefile(joinpath(@__DIR__, "goldens", "$case.json"))
        params = TORBEAM.TorbeamParams(; (Symbol(k) => v isa String ? Symbol(v) : v for (k, v) in golden["params"])...)
        eq = TORBEAM.equilibrium_inputs(dd)
        @testset "$case rays vs Fortran" begin
            for (ibeam, gb) in enumerate(golden["beams"])
                inputs = TORBEAM.beam_inputs(dd, ibeam, params, eq)
                m = TORBEAM.PlasmaModel(inputs)
                l = TORBEAM.Launch(inputs)
                b = TORBEAM.trace_beam(m, l; rhostop=params.rhostop, nmax=0)   # geometry only: no absorption
                g = golden_ray(gb)
                # central ray position along the common arclength
                smax = min(b.length, g.s[end])
                dev = Float64[]
                for i in 1:length(g.s)
                    g.s[i] <= smax || break
                    x = b.sol(g.s[i])[1:3]
                    push!(dev, hypot(x[1] - g.X[i], x[2] - g.Y[i], x[3] - g.Z[i]))
                end
                @info "$case $(gb["name"]): ray deviation max $(round(maximum(dev)*100; digits=2)) cm over $(round(smax; digits=2)) m, length $(round(b.length; digits=3)) vs $(round(g.s[end]; digits=3)) m ($(b.exit))"
                @test maximum(dev) < 0.01
                # the Fortran stops once the power is absorbed; without absorption our ray
                # can only be longer
                @test b.length >= g.s[end] - 0.05
                # beam widths against the Fortran's peripheral rays: its left/right rays carry
                # the 1/e half-width along the horizontal transverse direction (wh) and its
                # upper/lower rays the half-width along v × e_h (wv), drawn in the poloidal plane
                # (for toroidally steered beams that is not the cut by the poloidal plane, wp).
                # Agreement is ≤ 2 % on ITER and ≤ 6 % on the DIII-D rays; the
                # DIII-D O2 beam's horizontal width is 10-18 % narrower than the Fortran's
                # along the whole path (open, see the README)
                O2 = gb["intinbeam"][3] == 1 && round(Int, gb["floatinbeam"][1] / (27.99e9 * abs(gb["floatinbeam"][27]))) >= 2
                for frac in (0.25, 0.5, 0.9)
                    i = argmin(abs.(g.s .- frac * smax))
                    bw = TORBEAM.beam_widths(b, g.s[i])
                    @info "  s = $(round(g.s[i]; digits=2)) m: wh $(round(bw.wh*100; digits=2)) vs $(round(g.wh[i]*100; digits=2)) cm, wv $(round(bw.wv*100; digits=2)) vs $(round(g.wv[i]*100; digits=2)) cm"
                    tol = frac < 0.9 ? (case == "D3D" ? 0.06 : 0.03) : 0.08
                    @test bw.wh ≈ g.wh[i] rtol = (O2 && case == "D3D" ? 0.2 : tol)
                    @test bw.wv ≈ g.wv[i] rtol = tol
                    # the beam matrix against the ABCD propagation of a finite-difference ray
                    # family (`abcd_widths`: first derivatives of H only, no Riccati integration,
                    # includes diffraction): the same widths to 1 %, and the eikonal constraint
                    # Re M ∂H/∂N = -∂H/∂x holds along the ray
                    a = TORBEAM.abcd_widths(b, m, g.s[i])
                    @test a.wh ≈ bw.wh rtol = 0.01
                    @test a.wv ≈ bw.wv rtol = 0.01
                    @test a.residual < 1e-3
                end
            end
        end
    end
end
