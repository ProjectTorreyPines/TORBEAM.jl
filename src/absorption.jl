# Electron-cyclotron absorption along the ray (weak-damping approximation).
#
# The absorption coefficient of the power, α [1/m], follows from first-order
# perturbation of the cold dispersion relation by the anti-Hermitian part ε^a of
# the hot dielectric tensor: with Λ(N) = e*·D(N)·e, D = NN - N²I + ε_cold, and
# e the cold polarization (D e = 0), a damping δN = iκ v̂ along the ray gives
#     κ = -(e*·ε^a·e) / (v̂·∂Λ/∂N),   α = 2 k0 κ.
# ε^a is evaluated for a relativistic Maxwellian (Maxwell-Jüttner) exactly, as a
# line integral along the resonance curve γ = nY + N∥u∥ (an ellipse in the
# (u∥, u⊥) plane, u = p/mc) for harmonics n = 1..nmax.
#
# Local frame: x̂ along N⊥, ẑ along B, ŷ = ẑ × x̂. Electron charge sign enters
# through D = -XY/(1-Y²) in the cold tensor and the Bessel matrix below.

import SpecialFunctions: besselj, besselkx
import QuadGK: quadgk

"""
    cold_tensor(X, Y)

Stix parameters `(S, D, P)` of the cold electron plasma (`D < 0` for electrons)
"""
cold_tensor(X::Real, Y::Real) = (1 - X / (1 - Y^2), -X * Y / (1 - Y^2), 1 - X)

"""
    cold_polarization(X, Y, Nperp, Npar)

Unit polarization vector `e` (complex, local frame) of the cold wave with refractive
index `(Nperp, 0, Npar)`: the null vector of `D = NN - N²I + ε_cold`
"""
function cold_polarization(X::Real, Y::Real, Nperp::Real, Npar::Real)
    S, D, P = cold_tensor(X, Y)
    N2 = Nperp^2 + Npar^2
    M = [S-Npar^2 -im*D Nperp*Npar; im*D S-N2 0; Nperp*Npar 0 P-Nperp^2]
    best = zeros(ComplexF64, 3)
    bestn = 0.0
    for (i, j) in ((1, 2), (1, 3), (2, 3))
        c = cross(M[i, :], M[j, :])      # bilinear: orthogonal (without conjugation) to both rows
        nc = norm(c)
        if nc > bestn
            best = c / nc
            bestn = nc
        end
    end
    return best
end

"""
    bessel_matrix(n, b, upar, uperp)

Stix's harmonic matrix `T^n = w w†`, `w = (-(n/b)J_n, -iJ_n', J_n u∥/u⊥)` (electrons), for Bessel argument `b = N⊥ u⊥ / Y`
"""
function bessel_matrix(n::Int, b::Real, upar::Real, uperp::Real)
    Jn = besselj(n, b)
    Jnp = 0.5 * (besselj(n - 1, b) - besselj(n + 1, b))
    nJb = b > 1e-12 ? n * Jn / b : (n == 1 ? 0.5 : 0.0)     # (n/b) J_n(b), finite at b -> 0
    # T = w w† (positive semidefinite by construction). The sign of the first
    # component is the electron one: in the cold limit the n = ±1 terms must give
    # ε_xy = -iD with D = -XY/(1-Y²), which fixes it to -(n/b)J_n.
    w = [-nJb, -im * Jnp, Jn * upar / uperp]
    return w * w'
end

"""
    antihermitian_tensor(X, Y, Nperp, Npar, μ; nmax=3, rtol=1e-6)

Anti-Hermitian part of the relativistic Maxwellian dielectric tensor (local
frame), `μ = mc²/Te`, summed over harmonics `1:nmax`:
    ε^a = π X Σₙ ∫d³u (u⊥² μ f / γ) T^n δ(γ - nY - N∥u∥),  f = μ e^{-μγ} / (4π K₂(μ))
"""
function antihermitian_tensor(X::Real, Y::Real, Nperp::Real, Npar::Real, μ::Real; nmax::Int=3, rtol::Float64=1e-6)
    εa = zeros(ComplexF64, 3, 3)
    (X <= 0 || abs(Npar) >= 1) && return εa
    fnorm = μ / (4π * besselkx(2, μ))         # f = fnorm * exp(-μ(γ-1))  (besselkx = e^μ K₂(μ))
    den = 1 - Npar^2
    for n in 1:nmax
        A = n^2 * Y^2 + Npar^2 - 1              # resonance exists only for A > 0
        A > 0 || continue
        upc = n * Y * Npar / den
        apar = sqrt(A) / den
        aperp = sqrt(A / den)
        function integrand(φ)
            sφ, cφ = sincos(φ)
            upar = upc + apar * cφ
            uperp = aperp * sφ
            γ = n * Y + Npar * upar
            dl = hypot(apar * sφ, aperp * cφ)                      # curve element per dφ
            grad = hypot(uperp / γ, upar / γ - Npar)                 # |∇(γ - nY - N∥u∥)|
            w = 2π * uperp * uperp^2 * μ * fnorm * exp(-μ * (γ - 1)) / γ * dl / grad
            return w * bessel_matrix(n, Nperp * uperp / Y, upar, uperp)
        end
        I, _ = quadgk(integrand, 0.0, π; rtol, atol=0.0, order=11)
        εa .+= π * X .* I
    end
    return εa
end

"""
    absorption_coefficient(st, w::WaveParams, N, v; nmax=3)

Power absorption coefficient α [1/m] at a plasma `state` `st` for wave vector
`N` and ray direction `v` (both in the global Cartesian frame)
"""
function absorption_coefficient(st, w::WaveParams, N::AbstractVector, v::AbstractVector; nmax::Int=3)
    X, Y = plasma_XY(st, w)
    (X <= 0 || st.Te <= 0) && return 0.0
    bhat = collect(st.B) ./ st.Bmag
    Npar = dot(N, bhat)
    Nperp_vec = N .- Npar .* bhat
    Nperp = norm(Nperp_vec)
    if Nperp > 1e-10
        xhat = Nperp_vec ./ Nperp
    else
        xhat = cross(bhat, abs(bhat[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0])
        xhat ./= norm(xhat)
    end
    yhat = cross(bhat, xhat)
    μ = 510.99895 / st.Te                     # mc²/Te with Te in keV
    εa = antihermitian_tensor(X, Y, Nperp, Npar, μ; nmax)
    e = cold_polarization(X, Y, Nperp, Npar)
    num = real(dot(e, εa * e))                 # e*·ε^a·e  (dot conjugates its first argument)
    Nloc = [Nperp, 0.0, Npar]
    Ne = sum(Nloc .* e)                        # N·e, bilinear
    dΛ = [2 * real(conj(e[k]) * Ne) - 2 * Nloc[k] for k in 1:3]
    vhat = v ./ norm(v)
    vloc = [dot(vhat, xhat), dot(vhat, yhat), dot(vhat, bhat)]
    denom = dot(vloc, dΛ)
    κ = -num / denom
    return 2 * w.k0 * κ
end
