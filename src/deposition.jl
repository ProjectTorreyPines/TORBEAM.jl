# Power deposition profile: the power absorbed along the central ray is spread
# over the beam's Gaussian cross-section and binned in rho_pol; dP/dV follows
# from the flux-surface volumes computed on the equilibrium grid.

"""
    flux_volumes(m::PlasmaModel; nρ=201, refine=4)

Plasma volume [m³] inside the flux surfaces `ρ_pol = ρgrid` (uniform grid on
[0, 1]), by counting cells of the equilibrium grid refined `refine` times.
Only the region connected to the magnetic axis (through ψn < 0.995) counts, so
private-flux regions of a diverted equilibrium are excluded.
"""
function flux_volumes(m::PlasmaModel; nρ::Int=201, refine::Int=4)
    nR = refine * (length(m.R) - 1) + 1
    nZ = refine * (length(m.Z) - 1) + 1
    R = range(m.R[1], m.R[end]; length=nR)
    Z = range(m.Z[1], m.Z[end]; length=nZ)
    dA = step(R) * step(Z)
    ψn = [psi_norm(m, r, z) for r in R, z in Z]
    # connected component of the axis in ψn < 0.995, dilated by two cells
    core = ψn .< 0.995^2
    seen = falses(nR, nZ)
    i0 = clamp(round(Int, (m.R_axis - R[1]) / step(R)) + 1, 1, nR)
    j0 = clamp(round(Int, (m.Z_axis - Z[1]) / step(Z)) + 1, 1, nZ)
    stack = [(i0, j0)]
    seen[i0, j0] = true
    while !isempty(stack)
        i, j = pop!(stack)
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            if 1 <= ii <= nR && 1 <= jj <= nZ && !seen[ii, jj] && core[ii, jj]
                seen[ii, jj] = true
                push!(stack, (ii, jj))
            end
        end
    end
    inside = copy(seen)
    for _ in 1:2
        grown = copy(inside)
        for i in 1:nR, j in 1:nZ
            inside[i, j] && continue
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                ii, jj = i + di, j + dj
                if 1 <= ii <= nR && 1 <= jj <= nZ && inside[ii, jj]
                    grown[i, j] = true
                    break
                end
            end
        end
        inside = grown
    end
    ρgrid = range(0.0, 1.0; length=nρ)
    dV = zeros(nρ)   # volume in each ρ bin [ρ_k, ρ_{k+1})
    for i in 1:nR, j in 1:nZ
        (inside[i, j] && ψn[i, j] < 1) || continue
        ρ = sqrt(max(ψn[i, j], 0.0))
        k = clamp(floor(Int, ρ * (nρ - 1)) + 1, 1, nρ - 1)
        dV[k] += 2π * R[i] * dA
    end
    V = cumsum(dV)
    return ρgrid, [0.0; V[1:end-1]]
end

"""
    deposition(b::BeamSolution, m::PlasmaModel; nρ=NPNT, width_factor=1.0)

Absorbed power per unit volume `dP/dV` [W/m³] on the TORBEAM rho_pol grid
`range(0, 1, length=nρ+1)[1:nρ]`, together with the grid, the power absorbed
per bin and the volumes.

The power absorbed on each solver step is spread in rho_pol as a Gaussian whose
mean and width come from the beam cross-section: along each principal axis of
the power Gaussian (variance `σ² = 1/(2 k0 λ)`, `λ` eigenvalues of Φ restricted
to the plane ⟂ to the ray) rho_pol is sampled at ±σ, giving the linear
(gradient) and quadratic (curvature, e.g. near the axis or tangent to the flux
surfaces) contributions to its spread. Negative rho_pol is folded back.
"""
function deposition(b::BeamSolution, m::PlasmaModel; nρ::Int=NPNT, width_factor::Float64=1.0)
    ρgrid = range(0.0, 1.0; length=nρ + 1)[1:nρ]
    dρ = 1 / nρ
    Pbin = zeros(nρ)
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
        μρ = ρ0
        σρ2 = 0.0
        for a in 1:2
            # sampling offset capped at a quarter of the minor radius: beyond that the
            # linearized spread is meaningless anyway (beam wider than the plasma)
            σ = min(width_factor * sqrt(1 / (2 * k0 * max(ev.values[a], 1e-12))), 0.25 * m.a)
            d = ev.vectors[1, a] * e1 + ev.vectors[2, a] * e2
            xp = x0 + σ * d
            xm = x0 - σ * d
            ρp = rho_pol(m, hypot(xp[1], xp[2]), xp[3])
            ρm = rho_pol(m, hypot(xm[1], xm[2]), xm[3])
            σρ2 += (0.5 * (ρp - ρm))^2 + (0.5 * (ρp + ρm) - ρ0)^2
            μρ += 0.5 * (ρp + ρm) - ρ0
        end
        σρ = clamp(sqrt(σρ2), 0.5 * dρ, 0.5)
        μρ = clamp(μρ, 0.0, 1.5)
        # Gaussian in rho, folded at rho = 0, restricted to rho < 1
        klo = clamp(floor(Int, (μρ - 5σρ) / dρ), 1, nρ)
        khi = clamp(ceil(Int, (μρ + 5σρ) / dρ) + 1, 1, nρ)
        klo <= khi || continue
        wsum = 0.0
        w = zeros(khi - klo + 1)
        for (j, k) in enumerate(klo:khi)
            ρ = ρgrid[k] + 0.5 * dρ
            w[j] = exp(-0.5 * ((ρ - μρ) / σρ)^2) + exp(-0.5 * ((ρ + μρ) / σρ)^2)
            wsum += w[j]
        end
        wsum > 0 || continue
        for (j, k) in enumerate(klo:khi)
            Pbin[k] += ΔP * w[j] / wsum
        end
    end
    ρV, V = flux_volumes(m)
    Vc = cubic_resample(collect(ρV), V, collect(ρgrid))
    dVdρ = [k == 1 ? (Vc[2] - Vc[1]) / dρ : k == nρ ? (Vc[nρ] - Vc[nρ-1]) / dρ : (Vc[k+1] - Vc[k-1]) / (2dρ) for k in 1:nρ]
    dPdV = [dVdρ[k] > 0 ? Pbin[k] / dρ / dVdρ[k] : 0.0 for k in 1:nρ]
    return (; ρ=collect(ρgrid), dPdV, Pbin, V=Vc)
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
