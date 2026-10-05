# Power deposition profile: the power absorbed along the central ray is spread
# over the beam's Gaussian cross-section and binned in rho_pol; dP/dV follows
# from the flux-surface volumes computed on the equilibrium grid.

"""
    flux_volumes(m::PlasmaModel; nρ=101, nθ=360)

Plasma volume [m³] inside the flux surfaces `ρ_pol = ρgrid` (uniform grid on
[0, 1]): `V = π ∮ R² dZ` over the contour of each surface (`flux_contour`,
which follows the surfaces around the magnetic axis, so private-flux regions
of a diverted equilibrium are excluded). The separatrix itself is not traced;
`V(1)` is extrapolated quadratically in `ψn` from the last three surfaces.
"""
function flux_volumes(m::PlasmaModel; nρ::Int=101, nθ::Int=360)
    ρgrid = range(0.0, 1.0; length=nρ)
    V = zeros(nρ)
    for k in 2:nρ-1
        R, Z = flux_contour(m, ρgrid[k]; nθ)
        s = 0.0
        for i in 1:nθ
            ip = i == nθ ? 1 : i + 1
            s += 0.5 * (R[i]^2 + R[ip]^2) * (Z[ip] - Z[i])
        end
        V[k] = π * abs(s)
    end
    # V(ψn) is smooth through the edge: quadratic extrapolation in ψn = ρ²
    x = ρgrid[nρ-3:nρ-1] .^ 2
    y = V[nρ-3:nρ-1]
    x1 = 1.0
    V[nρ] = y[1] * (x1 - x[2]) * (x1 - x[3]) / ((x[1] - x[2]) * (x[1] - x[3])) +
            y[2] * (x1 - x[1]) * (x1 - x[3]) / ((x[2] - x[1]) * (x[2] - x[3])) +
            y[3] * (x1 - x[1]) * (x1 - x[2]) / ((x[3] - x[1]) * (x[3] - x[2]))
    return ρgrid, V
end

"""
    deposition(b::BeamSolution, m::PlasmaModel; nρ=NPNT, width_factor=1.0, efficiency=nothing, volumes=nothing, method=:maj, nsample=32)

Absorbed power per unit volume `dP/dV` [W/m³] and, when `efficiency(u, s)`
(local j∥/P_abs [A m/W] from the ray state) is given, the driven current
density `j` [A/m²] on the TORBEAM rho_pol grid `range(0, 1, length=nρ+1)[1:nρ]`,
together with the grid, the power and current×volume per bin and the volumes.

The power lost by the reference ray on each step is spread over the beam's
Gaussian amplitude (variances `σ² = 1/(2 k0 λ)` along the principal axes, `λ`
the eigenvalues of Φ restricted to the plane ⟂ to the ray; an `nsample`² grid over
±4.2σ) on a plane through the step's position, and each sample is assigned to
its own flux surface. The histogram is smoothed over a quarter of the step's
rho spread (estimated from samples at ±σ) with the kernel weighted by the bin
volumes, so the smoothing is in dP/dV and the bins next to the axis, whose
volume vanishes, do not blow up; dP/dV is the binned power over the exact bin
volume from `flux_volumes`. A part of the cross-section
displaced by ξ from the central ray is absorbed where its own path reaches the
resonance, not at the same arclength: for a beam crossing the resonance layer
obliquely this elongates the deposition region by 1/cos of the crossing angle.
`method=:maj` (the Fortran's `nprofcalc=1`, Poli et al. 2018 Eq. 14, whose
integration grid has one axis along the beam, one horizontal and one vertical):
the plane is the vertical one through the step, i.e. the resonance is taken as
a vertical surface, `δs = -(ξ·n)/(v̂·n)` with `n` the horizontal direction of the
ray. `method=:shifted` (`nprofcalc=2`): the plane is the local iso-Y surface,
`δs = -(ξ·∇Y)/(v̂·∇Y)`, limited to 2|ξ| for grazing crossings.
"""
function deposition(b::BeamSolution, m::PlasmaModel; nρ::Int=NPNT, width_factor::Float64=1.0, efficiency=nothing, volumes=nothing, method::Symbol=:maj, nsample::Int=32)
    method in (:maj, :shifted) || throw(ArgumentError("unknown deposition method $method"))
    ρgrid = range(0.0, 1.0; length=nρ + 1)[1:nρ]
    dρ = 1 / nρ
    Pbin = zeros(nρ)
    Jbin = zeros(nρ)      # driven current × dV, when `efficiency(u, s)` [A m/W] is given
    # volumes at the bin edges, interpolated in ψn = ρ² (V is linear in ψn near the axis)
    ρV, V = volumes === nothing ? flux_volumes(m) : volumes
    ρedges = collect(range(0.0, 1.0; length=nρ + 1))
    Ve = cubic_resample(collect(ρV) .^ 2, V, ρedges .^ 2)
    ΔV = max.(diff(Ve), 0.0)
    raw = zeros(nρ)
    work = zeros(nρ)
    tmp = zeros(nρ)
    k0 = b.launch.wave.k0
    ts = b.sol.t
    for i in 1:length(ts)-1
        s0, s1 = ts[i], ts[i+1]
        ΔP = power(b, s0) - power(b, s1)
        ΔP > 0 || continue
        s = 0.5 * (s0 + s1)
        u = b.sol(s)
        x0 = u[1:3]
        v = b.sol(s, Val{1})[1:3]
        v ./= norm(v)
        e1 = cross([0.0, 0.0, 1.0], v)
        e1 ./= norm(e1)
        e2 = cross(v, e1)
        Φ = imag(unpack_M(u))
        Φt = [dot(e1, Φ * e1) dot(e1, Φ * e2); dot(e2, Φ * e1) dot(e2, Φ * e2)]
        ev = eigen(Symmetric(Φt))
        ρ0 = rho_pol(m, hypot(x0[1], x0[2]), x0[3])
        # normal of the surface the samples are shifted onto: the gradient of the cyclotron
        # frequency (∝ |B|) for the iso-Y shift, the horizontal direction of the ray for the
        # vertical plane of the Fortran
        gY = method == :shifted ? ForwardDiff.gradient(p -> state(m, p[1], p[2], p[3]).Bmag, x0) : [v[1], v[2], 0.0]
        vY = dot(v, gY)
        shift(ξ) = abs(vY) > 1e-3 * norm(gY) ? -(dot(ξ, gY) / vY) : 0.0
        clampshift(δ, ξ) = method == :shifted ? clamp(δ, -2norm(ξ), 2norm(ξ)) : δ
        μρ = ρ0
        σρ2 = 0.0
        for a in 1:2
            # sampling offset capped at a quarter of the minor radius: beyond that the
            # linearized spread is meaningless anyway (beam wider than the plasma)
            σ = min(width_factor * sqrt(1 / (2 * k0 * max(ev.values[a], 1e-12))), 0.25 * m.a)
            d = ev.vectors[1, a] * e1 + ev.vectors[2, a] * e2
            δ = clampshift(shift(σ * d), σ * d)
            xp = x0 + σ * d + δ * v
            xm = x0 - σ * d - δ * v
            ρp = rho_pol(m, hypot(xp[1], xp[2]), xp[3])
            ρm = rho_pol(m, hypot(xm[1], xm[2]), xm[3])
            σρ2 += (0.5 * (ρp - ρm))^2 + (0.5 * (ρp + ρm) - ρ0)^2
            μρ += 0.5 * (ρp + ρm) - ρ0
        end
        σρ = clamp(sqrt(σρ2), 0.5 * dρ, 0.5)
        η = efficiency === nothing ? 0.0 : efficiency(u, s)
        # sample the power cross-section (uniform grid in the principal frame, Gaussian
        # weights), shift each sample along the ray to where its own path meets the
        # resonance, and bin the exact rho of each sample, smoothed over a fraction of the
        # step's rho spread (the kernel is weighted by the bin volumes: smoothing in dP/dV)
        tq, wq = sample_rule(nsample)
        d1 = ev.vectors[1, 1] * e1 + ev.vectors[2, 1] * e2
        d2 = ev.vectors[1, 2] * e1 + ev.vectors[2, 2] * e2
        σ1 = min(sqrt(1 / (2 * k0 * max(ev.values[1], 1e-12))), 0.5 * m.a)
        σ2 = min(sqrt(1 / (2 * k0 * max(ev.values[2], 1e-12))), 0.5 * m.a)
        kw = max(3.0, 0.25 * σρ / dρ)
        mbox = max(1, round(Int, kw))
        # raw histogram of the samples, each weighted by 1/D(k) with D = G ΔV the volume
        # under the (symmetric) smoothing operator G around its bin: P = ΔV ⊙ G r then
        # conserves the power exactly and smooths dP/dV rather than the binned power
        fill!(raw, 0.0)
        kmin, kmax = nρ, 1
        for (p, wp) in zip(tq, wq), (q, wqq) in zip(tq, wq)
            ξ = sqrt(2) * σ1 * p * d1 + sqrt(2) * σ2 * q * d2
            # (iso-Y: the shift is limited to twice the transverse displacement, beyond that
            # the straight continuation of a sub-ray is not reliable for grazing crossings)
            x = x0 + ξ + clampshift(shift(ξ), ξ) * v
            ρ = rho_pol(m, hypot(x[1], x[2]), x[3])
            ρ < 1 || continue
            kc = clamp(floor(Int, ρ / dρ) + 1, 1, nρ)
            raw[kc] += ΔP * wp * wqq
            kmin = min(kmin, kc)
            kmax = max(kmax, kc)
        end
        kmin <= kmax || continue
        lo = max(1, kmin - 6mbox)
        hi = min(nρ, kmax + 6mbox)
        # D = G ΔV on the window, then r = raw / D, then G r
        copyto!(work, lo, ΔV, lo, hi - lo + 1)
        box_smooth!(work, lo, hi, mbox, tmp)
        for k in lo:hi
            raw[k] = work[k] > 0 ? raw[k] / work[k] : 0.0
        end
        box_smooth!(raw, lo, hi, mbox, tmp)
        for k in lo:hi
            g = raw[k] * ΔV[k]
            Pbin[k] += g
            Jbin[k] += η * g
        end
    end
    dPdV = [ΔV[k] > 0 ? Pbin[k] / ΔV[k] : 0.0 for k in 1:nρ]
    j = [ΔV[k] > 0 ? Jbin[k] / ΔV[k] : 0.0 for k in 1:nρ]
    return (; ρ=collect(ρgrid), dPdV, j, Pbin, Jbin, V=Ve[1:nρ])
end

"""
    box_smooth!(y, lo, hi, m, tmp)

Three passes of a running mean of half-width `m` over `y[lo:hi]` in place
(truncated at the window ends without renormalisation): a Gaussian-like kernel of
standard deviation ≈ `m` bins that is symmetric as an operator, so
`∑ a ⊙ G b = ∑ b ⊙ G a` and the volume-weighted smoothing conserves the power exactly
"""
function box_smooth!(y::AbstractVector, lo::Int, hi::Int, m::Int, tmp::AbstractVector)
    for _ in 1:3
        # prefix sums on the window
        acc = 0.0
        for k in lo:hi
            acc += y[k]
            tmp[k] = acc
        end
        for k in lo:hi
            a = max(lo, k - m)
            b = min(hi, k + m)
            # fixed normalisation (mass beyond the window ends is dropped): keeps G symmetric
            y[k] = (tmp[b] - (a > lo ? tmp[a-1] : 0.0)) / (2m + 1)
        end
    end
    return y
end

"""
    sample_rule(n)

Transverse sampling of the power cross-section: `n` uniform nodes per axis over ±3 in
`t` (`ξ = √2 σ t`) with Gaussian weights normalised to one; a uniform grid bins more
smoothly than Gauss-Hermite nodes
"""
function sample_rule(n::Int)
    t = collect(range(-3.0, 3.0; length=n))
    w = exp.(-t .^ 2)
    return t, w ./ sum(w)
end

"""
    gauss_hermite(n)

Nodes and weights of the n-point Gauss-Hermite rule (weight exp(-t²)), from
the eigenvalues of the Jacobi matrix
"""
function gauss_hermite(n::Int)
    J = zeros(n, n)
    for i in 1:n-1
        J[i, i+1] = J[i+1, i] = sqrt(i / 2)
    end
    ev = eigen(Symmetric(J))
    return ev.values, sqrt(π) .* ev.vectors[1, :] .^ 2
end
