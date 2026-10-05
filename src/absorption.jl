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
import LinearAlgebra: dot, norm, cross, eigen, Symmetric, I, det

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
    harmonic_vector(n, b, upar, uperp)

Gyro-phase Fourier coefficient of the velocity, `w = ((n/b)J_n, iJ_n', J_n u∥/u⊥)`, for
Bessel argument `b = N⊥ u⊥ / Y`, in the local frame x̂ ∥ N⊥, ŷ = b̂ × x̂ with phase
e^{i(k·r - ωt)} and electrons gyrating counter-clockwise about b̂: along the orbit the
wave phase picks up `b (sin φ' - sin φ)` and `∫ (cos φ, sin φ, u∥/u⊥) e^{i(nφ - b sin φ)} dφ/2π`
gives `w`. The harmonic `n > 0` resonates at `γ = nY + N∥u∥`. In the cold limit the n = ±1
terms give ε_xy = -iD with D = -XY/(1-Y²); the x–z and y–z entries of `w w†` (odd in u∥)
are not visible there and follow from the same integral.
"""
function harmonic_vector(n::Int, b::Number, upar::Real, uperp::Real)
    Jn = besselj(n, b)
    Jnp = 0.5 * (besselj(n - 1, b) - besselj(n + 1, b))
    nJb = abs(b) > 1e-12 ? n * Jn / b : (abs(n) == 1 ? 0.5 * sign(n) : 0.0)     # (n/b) J_n(b), finite at b -> 0
    return [nJb, im * Jnp, Jn * upar / uperp]
end

"""
    bessel_matrix(n, b, upar, uperp)

Stix's harmonic matrix `T^n = w w†` (positive semidefinite by construction), see [`harmonic_vector`](@ref)
"""
function bessel_matrix(n::Int, b::Number, upar::Real, uperp::Real)
    return harmonic_matrix(harmonic_vector(n, b, upar, uperp))
end

"""
    harmonic_matrix(w)

`T = w w†` written as `w w̃ᵀ` with `w̃` the vector with the explicit `i` conjugated and the Bessel
functions not: identical for real Bessel argument, and the analytic continuation in `N⊥` (needed
for the complex root of the dispersion relation) otherwise
"""
harmonic_matrix(w::AbstractVector) = w * transpose([w[1], -w[2], w[3]])

"""
    antihermitian_tensor(X, Y, Nperp, Npar, μ; nmax=3, rtol=1e-6)

Anti-Hermitian part of the relativistic Maxwellian dielectric tensor (local
frame), `μ = mc²/Te`, summed over harmonics `1:nmax`:
    ε^a = π X Σₙ ∫d³u (u⊥² μ f / γ) T^n δ(γ - nY - N∥u∥),  f = μ e^{-μγ} / (4π K₂(μ))
"""
function antihermitian_tensor(X::Real, Y::Real, Nperp::Number, Npar::Real, μ::Real; nmax::Int=3, rtol::Float64=1e-6)
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
    absorption_coefficient(st, w::WaveParams, N, v; nmax=3, warm=true)

Power absorption coefficient α [1/m] at a plasma `state` `st` for wave vector
`N` and ray direction `v` (both in the global Cartesian frame). With `warm` the
perpendicular index, the polarization and the wave-matrix derivative come from
the warm dispersion relation where the plasma is resonant, and with `complex_root`
α = 2k₀ Im N⊥ (x̂·v̂) from the complex root of the full relativistic dispersion
relation (`complex_warm_root`, the Fortran's `nabsroutine=1`); otherwise the
weak-damping value α = 2k₀κ, κ = -(e*ε^a e)/(v̂·∂λ/∂N).
"""
function absorption_coefficient(st, w::WaveParams, N::AbstractVector, v::AbstractVector; nmax::Int=3, warm::Bool=true, complex_root::Bool=true)
    X, Y = plasma_XY(st, w)
    (X <= 0 || st.Te <= 0) && return 0.0
    ws = wave_state(st, w, N; nmax, warm)
    num = real(dot(ws.e, ws.εa * ws.e))            # e*·ε^a·e
    num > 0 || return 0.0
    # ∂λ/∂N in the local frame by central differences of λ(N⊥, N∥) (exact Hellmann-Feynman
    # for the cold case, and includes ∂ε^h/∂N for the warm one)
    h = 1e-4
    dλ_perp = (ws.λ(ws.Nperp + h, ws.Npar) - ws.λ(ws.Nperp - h, ws.Npar)) / (2h)
    dλ_par = (ws.λ(ws.Nperp, ws.Npar + h) - ws.λ(ws.Nperp, ws.Npar - h)) / (2h)
    vhat = v ./ norm(v)
    denom = dot(vhat, ws.xhat) * dλ_perp + dot(vhat, ws.bhat) * dλ_par
    denom < 0 || return 0.0                        # positive-energy wave: ∂λ/∂N antiparallel to v̂
    κ = -num / denom
    if warm && complex_root && ws.warm
        # TORBEAM's nabsroutine=1 route (Farina's WARMDISP): the imaginary part of N⊥ from the
        # complex root of the full relativistic dispersion relation, projected on the ray
        xv = dot(vhat, ws.xhat)
        if xv > 0.05
            z = complex_warm_root(ws.X, ws.Y, ws.Npar, ws.μ, w.mode, ws.Nperp, κ / xv; nmax)
            z === nothing || return 2 * w.k0 * imag(z) * xv
        end
    end
    return 2 * w.k0 * κ
end

# ---------------------------------------------------------------------------
# Warm dispersion: Hermitian part of the relativistic dielectric tensor, warm
# N⊥ and polarization (TORBEAM's nabsroutine = 1 route)
# ---------------------------------------------------------------------------
#
# ε^h = I + X Σₙ P∫ d³u (u⊥² f'/γ) Tⁿ / (γ - nY - N∥u∥). With E = nY + N∥u∥ the
# resonant denominator is rationalized, 1/(γ-E) = (γ+E)/Q, Q = γ² - E² a
# quadratic in u∥ whose real roots are the resonant u∥ (spurious zeros of γ+E for
# n < 0 carry a vanishing numerator). For each u⊥ the u∥ integral is split into
# a smooth part — the numerator minus its linear interpolant between the roots —
# done by Gauss-Legendre, and the principal value of the interpolant over Q,
# which is a closed form. The residues of the same roots give the anti-Hermitian
# part, so `hermitian_tensor` is checked against `antihermitian_tensor`.

"""
    gauss_legendre(n)

Nodes and weights of the n-point Gauss-Legendre rule on [-1, 1] (Golub-Welsch)
"""
function gauss_legendre(n::Int)
    J = zeros(n, n)
    for i in 1:n-1
        J[i, i+1] = J[i+1, i] = i / sqrt(4i^2 - 1)
    end
    ev = eigen(Symmetric(J))
    return ev.values, 2 .* ev.vectors[1, :] .^ 2
end

const GL_PAR = gauss_legendre(100)
const GL_PERP = gauss_legendre(32)

"""
    hermitian_tensor(X, Y, Nperp, Npar, μ; nmax=3)

Hermitian part (including the identity) of the relativistic Maxwellian
dielectric tensor in the local frame, harmonics `-1:nmax` (n ≤ -2 are O(b⁴))
"""
function hermitian_tensor(X::Real, Y::Real, Nperp::Number, Npar::Real, μ::Real; nmax::Int=3)
    εh = Matrix{ComplexF64}(I, 3, 3)
    (X <= 0 || abs(Npar) >= 1) && return εh
    fnorm = μ / (4π * besselkx(2, μ))
    c = 1 - Npar^2
    # velocity ranges where exp(-μ(γ-1)) matters
    γmax = 1 + 32 / μ
    L = min(sqrt(γmax^2 - 1), 3.0)
    tp, wp = GL_PAR
    tt, wt = GL_PERP
    for n in -1:nmax                                 # n ≤ -2 contribute O(b⁴) and are dropped
        A = n^2 * Y^2 + Npar^2 - 1
        Ures = A > 0 ? min(sqrt(A / c), L) : 0.0     # u⊥ range with real roots (resonance ellipse)
        acc = zeros(ComplexF64, 3, 3)
        # u⊥ nodes: [0, Ures] with u⊥ = Ures sin θ (clustered away from the ellipse tip),
        # then the non-resonant remainder [Ures, L] with plain Gauss-Legendre
        nodes = Tuple{Float64,Float64}[]
        if Ures > 0
            for (tθ, wθ) in zip(tt, wt)
                θ = (tθ + 1) * π / 4
                push!(nodes, (Ures * sin(θ), Ures * cos(θ) * π / 4 * wθ))
            end
        end
        if Ures < L
            for (t, wq) in zip(tt, wt)
                push!(nodes, (Ures + (L - Ures) * (t + 1) / 2, (L - Ures) / 2 * wq))
            end
        end
        for (uperp, dup) in nodes
            uperp > 1e-9 || continue
            b = Nperp * uperp / Y
            w1, w2, w3 = harmonic_vector(n, b, uperp, uperp)     # w3 = J_n (u∥/u⊥) at u∥ = u⊥
            # numerator N(u∥) = (u⊥² f'(γ)/γ)(γ+E) T, f' = -μ f_M, as a 3x3 Hermitian matrix
            function numer(upar)
                γ = sqrt(1 + uperp^2 + upar^2)
                E = n * Y + Npar * upar
                w = [w1, w2, w3 * upar / uperp]
                return (uperp^2 * (-μ * fnorm * exp(-μ * (γ - 1))) / γ * (γ + E)) .* harmonic_matrix(w)
            end
            Q(upar) = c * upar^2 - 2 * n * Y * Npar * upar + (1 + uperp^2 - n^2 * Y^2)
            disc = n^2 * Y^2 * Npar^2 - c * (1 + uperp^2 - n^2 * Y^2)
            Iu = zeros(ComplexF64, 3, 3)
            if disc > 0
                sd = sqrt(disc)
                u1 = (n * Y * Npar - sd) / c
                u2 = (n * Y * Npar + sd) / c
                N1 = numer(u1)
                N2 = numer(u2)
                for (t, wq) in zip(tp, wp)
                    upar = L * t
                    Ñ = N1 .+ (N2 .- N1) .* ((upar - u1) / (u2 - u1))
                    Iu .+= (L * wq) .* ((numer(upar) .- Ñ) ./ Q(upar))
                end
                # PV of the interpolant: partial fractions of Ñ / (c (u-u1)(u-u2))
                Iu .+= N1 ./ (c * (u1 - u2)) .* log(abs((L - u1) / (L + u1)))
                Iu .+= N2 ./ (c * (u2 - u1)) .* log(abs((L - u2) / (L + u2)))
            else
                for (t, wq) in zip(tp, wp)
                    upar = L * t
                    Iu .+= (L * wq) .* (numer(upar) ./ Q(upar))
                end
            end
            acc .+= (2π * uperp * dup) .* Iu
        end
        εh .+= X .* acc
    end
    return (εh .+ εh') ./ 2
end

"""
    warm_dispersion(X, Y, Npar, μ, mode, Nperp0; nmax=3)

Real `N⊥` solving `det(NN - N²I + ε^h(N⊥)) = 0` by Newton from the cold root
`Nperp0`; returns `(Nperp, Dw)` with the warm wave matrix, or the cold values if
the iteration fails
"""
function warm_dispersion(X::Real, Y::Real, Npar::Real, μ::Real, mode::Integer, Nperp0::Real; nmax::Int=3)
    function wavematrix(Np)
        N = [Np, 0.0, Npar]
        return N * N' .- (Np^2 + Npar^2) .* Matrix{ComplexF64}(I, 3, 3) .+ hermitian_tensor(X, Y, Np, Npar, μ; nmax)
    end
    f(Np) = real(det(wavematrix(Np)))
    # secant from the cold root
    Np0 = Nperp0
    Np1 = Nperp0 * 1.01 + 1e-4
    f0 = f(Np0)
    f1 = f(Np1)
    Np = Np1
    for _ in 1:10
        f1 == f0 && break
        δ = -f1 * (Np1 - Np0) / (f1 - f0)
        δ = clamp(δ, -0.2 * max(Np1, 0.1), 0.2 * max(Np1, 0.1))
        Np = max(Np1 + δ, 1e-6)
        Np0, f0 = Np1, f1
        Np1, f1 = Np, f(Np)
        abs(δ) < 1e-7 && break
    end
    if !(isfinite(Np) && 0 < Np < 3 && abs(Np - Nperp0) < 0.5 * max(Nperp0, 0.1))
        Np = Nperp0
    end
    return Np, wavematrix(Np)
end

"""
    complex_warm_root(X, Y, Npar, μ, mode, Nperp0, κ0; nmax=3)

Complex root `N⊥` of the full relativistic dispersion relation
`det(NN - N²I + ε^h(N⊥) + iε^a(N⊥)) = 0` at fixed real `N∥`, by a complex secant started at
`Nperp0 + iκ0` (the warm real root and the weak-damping imaginary part); the tensors are
continued analytically in `N⊥` through their Bessel arguments. Returns `nothing` when the
iteration does not converge to a nearby root with `Im N⊥ ≥ 0`.
"""
function complex_warm_root(X::Real, Y::Real, Npar::Real, μ::Real, mode::Integer, Nperp0::Real, κ0::Real; nmax::Int=3)
    function f(Np)
        N = [Np, 0.0, Npar]
        D = N * transpose(N) .- (Np^2 + Npar^2) .* Matrix{ComplexF64}(I, 3, 3) .+ hermitian_tensor(X, Y, Np, Npar, μ; nmax) .+ im .* antihermitian_tensor(X, Y, Np, Npar, μ; nmax)
        return det(D)
    end
    z0 = complex(Nperp0, κ0)
    z1 = complex(Nperp0, 1.5κ0 + 1e-5)
    f0, f1 = f(z0), f(z1)
    z = z1
    for _ in 1:12
        f1 == f0 && break
        δ = -f1 * (z1 - z0) / (f1 - f0)
        a = abs(δ)
        a > 0.1 * max(Nperp0, 0.1) && (δ *= 0.1 * max(Nperp0, 0.1) / a)
        z = z1 + δ
        z0, f0 = z1, f1
        z1, f1 = z, f(z)
        abs(δ) < 1e-8 && break
    end
    ok = isfinite(z) && imag(z) >= 0 && abs(real(z) - Nperp0) < 0.3 * max(Nperp0, 0.1) && abs(f1) < 1e-3 * abs(f(complex(Nperp0, κ0)))
    return ok ? z : nothing
end

"""
    null_vector(M)

Unit null vector of a (nearly) singular complex 3x3 matrix, from its row cross products
"""
function null_vector(M::AbstractMatrix)
    best = zeros(ComplexF64, 3)
    bestn = 0.0
    for (i, j) in ((1, 2), (1, 3), (2, 3))
        c = cross(M[i, :], M[j, :])
        nc = norm(c)
        if nc > bestn
            best = c / nc
            bestn = nc
        end
    end
    return best
end

"""
    wave_state(st, w::WaveParams, N; nmax=3, warm=true)

Local wave quantities at a plasma `state`: the frame (`xhat`, `yhat`, `bhat`),
`Nperp`, `Npar`, the polarization `e`, the anti-Hermitian tensor `εa` and the
function `λ(N⊥, N∥) = e*·D·e` used for the group-velocity derivative. With `warm`
the perpendicular index and the polarization come from the warm dispersion
relation (Hermitian relativistic tensor) wherever the plasma is resonant.
"""
function wave_state(st, w::WaveParams, N::AbstractVector; nmax::Int=3, warm::Bool=true)
    X, Y = plasma_XY(st, w)
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
    μ = 510.99895 / max(st.Te, 1e-3)
    e = cold_polarization(X, Y, Nperp, Npar)
    εa = antihermitian_tensor(X, Y, Nperp, Npar, μ; nmax)
    S, D, P = cold_tensor(X, Y)
    εcold = [S -im*D 0; im*D S 0; 0 0 P]
    λ(Np, Npl) = (Nv = [Np, 0.0, Npl]; real(dot(e, (Nv * Nv' - (Np^2 + Npl^2) * I + εcold) * e)))
    resonant = real(dot(e, εa * e)) > 1e-8
    if warm && resonant
        Nperp_w, Dw = warm_dispersion(X, Y, Npar, μ, w.mode, Nperp; nmax)
        e = null_vector(Dw)
        εa = antihermitian_tensor(X, Y, Nperp_w, Npar, μ; nmax)
        Nperp = Nperp_w
        λ = (Np, Npl) -> (Nv = [Np, 0.0, Npl]; real(dot(e, (Nv * Nv' - (Np^2 + Npl^2) * I + hermitian_tensor(X, Y, Np, Npl, μ; nmax)) * e)))
    end
    return (; xhat, yhat, bhat, Nperp, Npar, e, εa, λ, X, Y, μ, warm=warm && resonant)
end
