# Electron-cyclotron current drive by the adjoint (Spitzer-function) method.
#
# The driven parallel current per absorbed power is
#     j∥ / P_abs = -(e / m_e c) ⟨ d·∇_u χ ⟩_W ,
# the average over the resonance curve, with the same weights W as the
# absorption, of the derivative of the response function χ along the
# quasilinear diffusion direction d = (N∥, (γ - N∥u∥)/u⊥) in (u∥, u⊥).
#
# χ solves the bounce-averaged adjoint equation with the high-velocity
# relativistic test-particle collision operator (slowing down on electrons,
# pitch-angle scattering on electrons and ions with Z_eff):
#     -ν_s u ∂χ/∂u + (ν_d/2) ⟨L_ξ⟩_b χ = -(u/γ) ⟨ξ⟩_b    (source: the parallel velocity),
#     ν_s = ν₀ γ²/u³,  ν_d = ν₀ (1 + Z_eff) γ/u³,  ν₀ = n_e e⁴ lnΛ / (4π ε₀² m_e² c³).
# ν_s is the drag on the speed (the friction on the velocity vector is twice
# that, but the perpendicular diffusion returns half of the energy) and
# (ν_d/2) L = ⟨Δv⊥²⟩/(4v²) L is the pitch-angle part of the Fokker-Planck
# operator. In a uniform field this gives the classic χ = u⁴ξ / (ν₀(5+Z)). On a flux
# surface, with λ = u⊥²/(u² b), b = B/B_min, the pitch-angle operator becomes
# 4 ∂_λ[λ I(λ) ∂_λ] / J(λ) with the bounce integrals I = ∮ dl |ξ|/b and
# J = ∮ dl/|ξ| along the field line, χ vanishes for trapped electrons
# (λ > λ_c = B_min/B_max), and the equation is marched in u from χ(0) = 0 with
# an implicit tridiagonal solve in λ. χ̂ = ν₀ χ is stored (dimensionless).
#
# Besides this exact 2-D solution (`SpitzerFunction(fs, Zeff)`, `ncdroutine=3`)
# and the full linearized operator (`full_operator_response`, `ncdroutine=4`),
# the separable model of Lin-Liu, Chan & Prater, Phys. Plasmas 10 (2003) 4064,
# that the Fortran uses (`linliu_response`, `ncdroutine=1`) is available: the
# slowing-down term is kept to its l = 1 Legendre moment, which makes
# χ̂ = sgn(u∥) F(u) H(λ) separable with the circulating fraction f_c in F.

import SparseArrays

"""
    FluxSurface(m::PlasmaModel, ρ; nθ=180)

Flux surface ρ_pol = ρ traced around the magnetic axis: `R`, `Z`, `B` [T],
field-line length weights `dl` (= B/B_p dl_pol), `Bmin`, `Bmax`, trapped-passing
boundary `λc = Bmin/Bmax`
"""
struct FluxSurface
    ρ::Float64
    R::Vector{Float64}
    Z::Vector{Float64}
    B::Vector{Float64}
    dl::Vector{Float64}
    Bmin::Float64
    Bmax::Float64
    λc::Float64
end

function FluxSurface(m::PlasmaModel, ρ::Real; nθ::Int=180)
    ψn_target = ρ^2
    θs = range(0, 2π; length=nθ + 1)[1:nθ]
    R = zeros(nθ)
    Z = zeros(nθ)
    for (i, θ) in enumerate(θs)
        cθ, sθ = cos(θ), sin(θ)
        f(r) = psi_norm(m, m.R_axis + r * cθ, m.Z_axis + r * sθ) - ψn_target
        # search up to the edge of the equilibrium grid along this direction
        rmax = min(cθ > 0 ? (m.R[end] - m.R_axis) / cθ : cθ < 0 ? (m.R[1] - m.R_axis) / cθ : Inf,
                   sθ > 0 ? (m.Z[end] - m.Z_axis) / sθ : sθ < 0 ? (m.Z[1] - m.Z_axis) / sθ : Inf) * 0.999
        # bracket the first crossing from the axis outwards, then bisect
        lo, hi = 0.0, rmax
        n = 400
        for k in 1:n
            r = rmax * k / n
            if f(r) > 0
                lo, hi = rmax * (k - 1) / n, r
                break
            end
        end
        for _ in 1:60
            mid = 0.5 * (lo + hi)
            f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        r = 0.5 * (lo + hi)
        R[i] = m.R_axis + r * cθ
        Z[i] = m.Z_axis + r * sθ
    end
    B = zeros(nθ)
    Bp = zeros(nθ)
    for i in 1:nθ
        BR, Bφ, BZ = B_cyl(m, R[i], Z[i])
        B[i] = sqrt(BR^2 + Bφ^2 + BZ^2)
        Bp[i] = hypot(BR, BZ)
    end
    # poloidal arc elements (centred) times B/B_p
    dl = zeros(nθ)
    for i in 1:nθ
        ip = i == nθ ? 1 : i + 1
        im_ = i == 1 ? nθ : i - 1
        dl[i] = 0.5 * hypot(R[ip] - R[im_], Z[ip] - Z[im_]) * B[i] / max(Bp[i], 1e-6)
    end
    Bmin, Bmax = minimum(B), maximum(B)
    return FluxSurface(ρ, R, Z, B, dl, Bmin, Bmax, Bmin / Bmax)
end

"""
    bounce_integrals(fs::FluxSurface, λ)

`I(λ) = ∮ dl |ξ|/b`, `J(λ) = ∮ dl/|ξ|` and `L = ∮ dl` for `1 - ξ² = λ b`, `b = B/Bmin`
"""
function bounce_integrals(fs::FluxSurface, λ::Real)
    I = 0.0
    J = 0.0
    L = 0.0
    for i in eachindex(fs.B)
        b = fs.B[i] / fs.Bmin
        ξ2 = 1 - λ * b
        ξ2 > 0 || continue
        ξ = sqrt(ξ2)
        I += fs.dl[i] * ξ / b
        J += fs.dl[i] / ξ
        L += fs.dl[i]
    end
    return I, J, L
end

"""
    SpitzerFunction

`χ̂(u, λ) = ν₀ χ` on a flux surface (odd in ξ: the stored function is for ξ > 0),
as a cubic spline, plus the surface
"""
struct SpitzerFunction{S}
    fs::FluxSurface
    χ::S
    umax::Float64
    scale::Float64   # multiplies the efficiency: ⟨B/B_max⟩ for the Lin-Liu model (its ⟨j∥⟩ definition), 1 otherwise
end
SpitzerFunction(fs::FluxSurface, χ, umax::Real) = SpitzerFunction(fs, χ, Float64(umax), 1.0)

"""
    SpitzerFunction(fs::FluxSurface, Zeff; nu=300, nλ=200, umax=1.5)

Lorentz-model response: solve the bounce-averaged adjoint equation (see the
file header) on the surface by marching in u
"""
function SpitzerFunction(fs::FluxSurface, Zeff::Real; nu::Int=300, nλ::Int=200, umax::Float64=1.5)
    us = range(0.0, umax; length=nu)
    # λ grid clustered towards the trapped boundary, where χ ∝ sqrt(λc - λ)
    ts = range(0.0, 1.0; length=nλ)
    λs = fs.λc .* (1 .- (1 .- ts) .^ 2)
    λh = [0.5 * (λs[k] + λs[k+1]) for k in 1:nλ-1]         # half points
    Ih = [bounce_integrals(fs, λ)[1] for λ in λh]
    JL = [bounce_integrals(fs, λ) for λ in λs[1:nλ-1]]
    J = [x[2] for x in JL]
    L = [x[3] for x in JL]
    ξbar = L ./ J                                   # bounce-averaged |ξ| (passing)
    χ = zeros(nu, nλ)                              # χ(:, nλ) = 0: trapped boundary
    # operator (4/J_k) ∂_λ[λ I ∂_λ χ] on the interior nodes with nonuniform
    # differences; node 1 (λ = 0): no flux through λ = 0
    n = nλ - 1                                     # unknowns: nodes 1..nλ-1
    am = zeros(n)
    ap = zeros(n)
    for k in 1:n
        hp = λs[k+1] - λs[k]
        hm = k == 1 ? hp : λs[k] - λs[k-1]
        hc = 0.5 * (hp + hm)
        ap[k] = λh[k] * Ih[k] / (hp * hc)
        am[k] = k == 1 ? 0.0 : λh[k-1] * Ih[k-1] / (hm * hc)
    end
    for iu in 2:nu
        u = us[iu]
        du = us[iu] - us[iu-1]
        γ = sqrt(1 + u^2)
        νs = γ^2 / u^3
        νd = (1 + Zeff) * γ / u^3
        # implicit step: (χ_new - χ_old)/du = [ (νd/2) A χ_new + (u/γ) ξbar ] / (νs u); the adjoint
        # source is the parallel velocity u∥/γ (the current is -e∫v∥f₁), not the momentum
        c = du / (νs * u)
        lower = zeros(n - 1)
        diag = zeros(n)
        upper = zeros(n - 1)
        rhs = zeros(n)
        for k in 1:n
            a = 4 / J[k] * (νd / 2)
            diag[k] = 1 + c * a * (am[k] + ap[k])
            k > 1 && (lower[k-1] = -c * a * am[k])
            k < n && (upper[k] = -c * a * ap[k])        # χ at node nλ is 0
            rhs[k] = χ[iu-1, k] + c * u / γ * ξbar[k]
        end
        χ[iu, 1:n] = tridiagonal_solve(lower, diag, upper, rhs)
    end
    # spline on the (u, t) grid; evaluation maps λ -> t
    spl = Interpolations.cubic_spline_interpolation((us, ts), χ; extrapolation_bc=Interpolations.Line())
    return SpitzerFunction(fs, spl, umax)
end

"""
    flux_average(fs::FluxSurface, A)

Flux-surface average `∮ A dl_p/B_p / ∮ dl_p/B_p` of the samples `A` on the surface points
"""
flux_average(fs::FluxSurface, A::AbstractVector) = sum(A .* fs.dl ./ fs.B) / sum(fs.dl ./ fs.B)

"""
    circulating_fraction(fs::FluxSurface)

Effective circulating fraction `f_c = ¾ ⟨B²/B²max⟩ ∫₀¹ λ dλ / ⟨(1 - λ B/Bmax)^{1/2}⟩`
(Lin-Liu et al. 2003 Eq. 32; `1 - f_c` is the effective trapped fraction)
"""
function circulating_fraction(fs::FluxSurface)
    b = fs.B ./ fs.Bmax
    # λ = 1 - s²: removes the integrable 1/√(1-λ) singularity of the uniform-B limit
    ns = 2000
    g = [(sλ = (k - 0.5) / ns; λ = 1 - sλ^2; 2sλ * λ / flux_average(fs, sqrt.(max.(1 .- λ .* b, 0.0)))) for k in 1:ns]
    return 0.75 * flux_average(fs, b .^ 2) * sum(g) / ns
end

"""
    linliu_response(fs::FluxSurface, Zeff, μ; nu=300, nλ=200, umax=1.5)

The separable response of Lin-Liu, Chan & Prater (2003), Eqs. 27-33, as a `SpitzerFunction`
in the normalisation of this file (`χ̂ = u_e⁴ χ̃` so that the uniform non-relativistic limit is
`u⁴ξ/(5+Z)`): `χ̂ = sgn(u∥) F̂(u) H(λ)` with
`H(λ) = ½ ∫_λ^1 dλ'/⟨(1 - λ' B/Bmax)^{1/2}⟩` for `λ = (Bmax/B) u⊥²/u² ≤ 1` (0 for trapped) and
`F̂(u) = (u⁴/f_c) ∫₀¹ x^{ρ̂+3} (1 + u²x²)^{-3/2} [(1 + γ)/(1 + γ(ux))]^{ρ̂} dx`, `ρ̂ = (Z+1)/f_c`.
The slowing-down term is kept to its l = 1 Legendre moment (exact in the Lorentz-gas limit);
`scale = ⟨B/Bmax⟩` carries the prefactor of Lin-Liu's efficiency (Eq. 38), whose current is the
flux-surface average `⟨j∥⟩`.
"""
function linliu_response(fs::FluxSurface, Zeff::Real, μ::Real; nu::Int=300, nλ::Int=200, umax::Float64=1.5, momentum_conservation::Bool=false)
    fc = circulating_fraction(fs)
    ρ̂ = (Zeff + 1) / fc
    bmax = fs.B ./ fs.Bmax
    # H on the file's t grid: λ_code = λc (1 - (1-t)²) ⇔ λ_LL = λ_code/λc = 1 - (1-t)², and with
    # λ' = 1 - s'² the integral is H = ∫₀^{1-t} s' ds' / ⟨(1 - (1-s'²) B/Bmax)^{1/2}⟩ (no singularity)
    ts = range(0.0, 1.0; length=nλ)
    nfine = 4000
    sfine = range(0.0, 1.0; length=nfine + 1)
    gs = [(sλ = (k - 0.5) / nfine; sλ / flux_average(fs, sqrt.(max.(1 .- (1 - sλ^2) .* bmax, 0.0)))) for k in 1:nfine]
    cum = vcat(0.0, cumsum(gs) ./ nfine)            # ∫₀^{s} g ds' at the nodes sfine
    H = [IMAS.interp1d(collect(sfine), cum).(1 - t) for t in ts]
    us = range(0.0, umax; length=nu)
    # F̂ of the high-speed limit on the u grid (Eq. 33), Gauss-Legendre in x
    xg, wg = gauss_legendre(64)
    function Fhsl(u)
        u > 0 || return 0.0
        γ = sqrt(1 + u^2)
        acc = 0.0
        for (x, w) in zip(xg, wg)
            xx = 0.5 * (x + 1)
            γx = sqrt(1 + (u * xx)^2)
            acc += 0.5 * w * xx^(ρ̂ + 3) * (1 + (u * xx)^2)^(-1.5) * ((1 + γ) / (1 + γx))^ρ̂
        end
        return u^4 / fc * acc
    end
    F = Fhsl.(us)
    if momentum_conservation
        # Momentum conservation as the ratio R(x) = K_mc(x)/K_hsl(x) of the non-relativistic
        # l = 1 solutions with the trapped-particle sink (Romé et al. 1998, Appendix: the variational
        # polynomial `variational_spitzer` for x = u/u_e ≤ 3, where it is within 3 % of the exact 1-D
        # solution `SpitzerFunction1D(Z; sink)`, which takes over beyond and gives R → 1), applied to the fully
        # relativistic high-speed-limit F̂ of Eq. 33 — the "relativistic adaptation" of that
        # Appendix. (Marushchenko's weakly relativistic μ⁻¹ expansion of the polynomial, also in
        # `variational_spitzer`, is not used here: against the fully relativistic high-speed limit
        # its ratio grows with x instead of tending to 1, and the Fortran's enhancement does not.)
        uT = sqrt(2 / μ)
        _, χa = variational_spitzer(fc, Zeff, Inf)
        sp = SpitzerFunction1D(Zeff; sink=(1 - fc) / fc)
        xmax = 3.0
        for (i, u) in enumerate(us)
            x = u / uT
            Kmc = x <= xmax ? χa(x) : spitzer_ratio_sink(sp, x)          # non-relativistic K_mc(x)
            # F̂_hsl / K_hsl,nr = (u_e⁴/f_c) × relativistic factor of Eq. 33 (→ 1 as u → 0)
            rel = u > 0 ? F[i] / (u^4 / (fc * (4 + ρ̂))) : 1.0
            F[i] = uT^4 / fc * rel * Kmc
        end
    end
    χ = F * H'
    spl = Interpolations.cubic_spline_interpolation((us, ts), χ; extrapolation_bc=Interpolations.Line())
    return SpitzerFunction(fs, spl, umax, flux_average(fs, bmax))
end

"""
    full_operator_response(fs::FluxSurface, Zeff, μ; nu=200, nλ=150, umax=1.5, nbasis=30, field=true)

Full-operator response (a `SpitzerFunction`) on the surface, `μ = mc²/Te` (see the file header).
`nbasis` hat functions represent the l = 1 moment in the Schur complement;
`field=false` keeps only the test-particle part.
"""
function full_operator_response(fs::FluxSurface, Zeff::Real, μ::Real; nu::Int=200, nλ::Int=150, umax::Float64=1.5, nbasis::Int=30, field::Bool=true)
    uT = sqrt(2 / μ)                      # thermal velocity v_T/c
    c3 = 1 / uT^3                         # ν̂/ν₀
    us = collect(range(0.0, umax; length=nu + 1))[2:end]
    du = us[2] - us[1]
    ts = range(0.0, 1.0; length=nλ)
    λs = fs.λc .* (1 .- (1 .- ts) .^ 2)
    λh = [0.5 * (λs[k] + λs[k+1]) for k in 1:nλ-1]
    Ih = [bounce_integrals(fs, λ)[1] for λ in λh]
    JL = [bounce_integrals(fs, λ) for λ in λs[1:nλ-1]]
    J = [x[2] for x in JL]
    L = [x[3] for x in JL]
    ξbar = L ./ J
    nk = nλ - 1
    am = zeros(nk)
    ap = zeros(nk)
    for k in 1:nk
        hp = λs[k+1] - λs[k]
        hm = k == 1 ? hp : λs[k] - λs[k-1]
        hc = 0.5 * (hp + hm)
        ap[k] = λh[k] * Ih[k] / (hp * hc)
        am[k] = k == 1 ? 0.0 : λh[k-1] * Ih[k-1] / (hm * hc)
    end
    # exact thermal rates (Chandrasekhar) with the relativistic γ factors of the high-velocity
    # limit: slowing-down νs = c3 (2G/x) γ², pitch-angle νD = c3 γ (φ - G + Z)/x³ (as in the Lorentz
    # march), and the energy diffusion fixed by detailed balance on the relativistic Maxwellian,
    # D = ½νpar u² = νs γ/μ, i.e. νpar = c3 (2G/x³) γ³ — with the γ of the non-relativistic form
    # the drag and the diffusion acting on f_M no longer cancel and the leftover νs(1 - 1/γ²)
    # acts as a spurious drag on suprathermal electrons
    function rates(u)
        x = u / uT
        γ = sqrt(1 + u^2)
        φ = erf(x)
        G = chandrasekhar(x)
        return (νpar=c3 * 2G / x^3 * γ^3, νD=c3 * γ * (φ - G + Zeff) / x^3)
    end
    idx(i, k) = (k - 1) * nu + i
    N = nu * nk
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for k in 1:nk, i in 1:nu
        r = idx(i, k)
        u = us[i]
        γ = sqrt(1 + u^2)
        # u part of C(f_M χ)/f_M: (u³P)'/u² - (μu²/γ)P with P = ½νpar u χ' — the drag νs χ and the
        # diffusion acting on f_M, -½νpar u (μu/γ) χ, cancel identically by detailed balance, so the
        # operator is pure diffusion; conservative on the half points, χ₀ = 0 at u = 0, ghost χ ∝ u⁴
        # beyond umax
        coef = zeros(3)                  # offsets -1, 0, +1
        for sgn in (+1, -1)
            uh = u + sgn * du / 2
            rh = rates(uh)
            a_i = -0.5 * rh.νpar * uh * sgn / du
            a_n = 0.5 * rh.νpar * uh * sgn / du
            w = sgn * uh^3 / (du * u^2) - (μ * u^2 / γ) / 2
            coef[2] += w * a_i
            coef[2+sgn] += w * a_n
        end
        for (off, v) in zip((-1, 0, 1), coef)
            j = i + off
            if j == 0
                continue
            elseif j == nu + 1
                push!(rows, r); push!(cols, idx(i, k)); push!(vals, v * (1 + 4du / u))
            else
                push!(rows, r); push!(cols, idx(j, k)); push!(vals, v)
            end
        end
        # pitch-angle part: (νD/2)(4/J_k)[am χ_{k-1} - (am+ap) χ_k + ap χ_{k+1}]
        a = 4 / J[k] * rates(u).νD / 2
        push!(rows, r); push!(cols, r); push!(vals, -a * (am[k] + ap[k]))
        k > 1 && (push!(rows, r); push!(cols, idx(i, k - 1)); push!(vals, a * am[k]))
        k < nk && (push!(rows, r); push!(cols, idx(i, k + 1)); push!(vals, a * ap[k]))
    end
    A = SparseArrays.sparse(rows, cols, vals, N, N)
    Alu = SparseArrays.lu(A)
    source(f::AbstractVector) = vec([-ξbar[k] * f[i] for i in 1:nu, k in 1:nk])
    χv = Alu \ source(us ./ sqrt.(1 .+ us .^ 2))     # adjoint source: parallel velocity u∥/γ
    if field
        # l = 1 moment at the outboard point: g₁(u) = (3/2)∫₀^λc χ dλ
        wλ = zeros(nk)
        for k in 1:nk
            lo = k == 1 ? λs[1] : 0.5 * (λs[k-1] + λs[k])
            hi = 0.5 * (λs[k] + λs[k+1])
            wλ[k] = 1.5 * (hi - lo)
        end
        moment(v) = reshape(v, nu, nk) * wλ
        # with χ constant along the orbit the local moment is g₁(θ) = b(θ) g₁(outboard) (the
        # parallel flow follows B), so the field term ξ(θ) F₁g₁(θ) bounce-averages to
        # ⟨ξ b⟩ F₁ g₁(outboard) = c_b ξ̄ F₁ g₁(outboard) with c_b = ∮ b dl / ∮ dl
        cb = sum(fs.B ./ fs.Bmin .* fs.dl) / sum(fs.dl)
        F1 = (cb * c3) .* field_matrix(us ./ uT)
        # Schur complement on nb hat functions: g₁ = Σ c_j B_j, χ = χ₀ + Σ c_j Ψ_j,
        # Ψ_j = A⁻¹ source(F₁ B_j), c from collocation of g₁ = moment(χ) at the coarse nodes
        nb = min(nbasis, nu)
        cj = round.(Int, range(1, nu; length=nb))
        B = zeros(nu, nb)
        for j in 1:nb
            lo = j == 1 ? cj[1] : cj[j-1]
            hi = j == nb ? cj[nb] : cj[j+1]
            for i in lo:hi
                B[i, j] = i <= cj[j] ? (lo == cj[j] ? 1.0 : (i - lo) / (cj[j] - lo)) : (hi - i) / (hi - cj[j])
            end
        end
        Ψ = zeros(N, nb)
        for j in 1:nb
            Ψ[:, j] = Alu \ source(F1 * B[:, j])
        end
        Kc = zeros(nb, nb)
        for j in 1:nb
            Kc[:, j] = moment(Ψ[:, j])[cj]
        end
        c = (I - Kc) \ moment(χv)[cj]
        χv += Ψ * c
    end
    χ = zeros(nu + 1, nλ)                 # u = 0 row and λ = λc column are 0
    χ[2:end, 1:nk] = reshape(χv, nu, nk)
    spl = Interpolations.cubic_spline_interpolation((range(0.0, umax; length=nu + 1), ts), χ; extrapolation_bc=Interpolations.Line())
    return SpitzerFunction(fs, spl, umax)
end

"""
    tridiagonal_solve(lower, diag, upper, rhs)

Thomas algorithm
"""
function tridiagonal_solve(lower::AbstractVector, diag::AbstractVector, upper::AbstractVector, rhs::AbstractVector)
    n = length(diag)
    c = copy(upper)
    d = copy(rhs)
    b = copy(diag)
    for i in 2:n
        w = lower[i-1] / b[i-1]
        b[i] -= w * c[i-1]
        d[i] -= w * d[i-1]
    end
    x = zeros(n)
    x[n] = d[n] / b[n]
    for i in n-1:-1:1
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    end
    return x
end

"""
    chi(sf::SpitzerFunction, u, λ)

`χ̂(u, λ)` (for ξ > 0)
"""
chi(sf::SpitzerFunction, u::Real, λ::Real) = sf.χ(min(u, sf.umax), 1 - sqrt(max(1 - min(λ, sf.fs.λc) / sf.fs.λc, 0.0)))

"""
    chi_derivatives(sf::SpitzerFunction, upar, uperp, b)

`d·∇_u χ̂` for the diffusion direction `d` is assembled by the caller from
`(∂χ̂/∂u, ∂χ̂/∂λ)` returned here, at `(u∥, u⊥)` on a point of the surface where
`b = B/Bmin`; `χ̂` is odd in `u∥`
"""
function chi_derivatives(sf::SpitzerFunction, upar::Real, uperp::Real, b::Real)
    u = hypot(upar, uperp)
    λc = sf.fs.λc
    λ = min(uperp^2 / (u^2 * b), λc)
    t = 1 - sqrt(max(1 - λ / λc, 0.0))               # grid variable: λ = λc (1 - (1-t)²)
    g = Interpolations.gradient(sf.χ, min(u, sf.umax), t)
    dtdλ = t < 1 ? 1 / (2λc * (1 - t)) : 0.0
    σ = sign(upar)
    return σ * g[1], σ * g[2] * dtdλ, λ
end

"""
    cd_efficiency(sf::SpitzerFunction, b, st, w::WaveParams, N; nmax=3, warm=true)

Local driven current per absorbed power, `j∥/P_abs` [A/W per m² → A m / W],
positive along B, at a point with `b = B/Bmin` on the surface of `sf`:
`-(e/(m_e c ν₀)) ⟨d·∇_u χ̂⟩_W` with ν₀ from the local density and the Coulomb
logarithm; the weights use the (warm, if `warm`) polarization and N⊥
"""
function cd_efficiency(sf::SpitzerFunction, b::Real, st, w::WaveParams, N::AbstractVector; nmax::Int=3, warm::Bool=true)
    X, Y = plasma_XY(st, w)
    (X <= 0 || st.Te <= 0) && return 0.0
    ws = wave_state(st, w, N; nmax, warm)
    num, den = cd_resonance(X, Y, ws.Nperp, ws.Npar, ws.μ, ws.e, sf, b; nmax)
    den > 0 || return 0.0
    lnΛ = coulomb_log(st.ne, st.Te)
    ν0 = st.ne * e_charge^4 * lnΛ / (4π * ε_0^2 * m_e^2 * c_light^3)
    return -sf.scale * e_charge / (m_e * c_light * ν0) * num / den
end

"""
    coulomb_log(ne, Te)

Electron Coulomb logarithm (NRL), `ne` [m^-3], `Te` [keV]
"""
coulomb_log(ne::Real, Te::Real) = 24 - log(sqrt(ne * 1e-6) / (Te * 1e3))

"""
    cd_resonance(X, Y, Nperp, Npar, μ, e, sf, b; nmax, rtol)

Resonance-curve integrals `∫ W (d·∇χ̂)` and `∫ W` with the absorption weight
`W ∝ (u⊥² μ f/γ) |w†e|² δ(γ - nY - N∥u∥)`
"""
function cd_resonance(X::Real, Y::Real, Nperp::Real, Npar::Real, μ::Real, e::AbstractVector, sf::SpitzerFunction, b::Real; nmax::Int=3, rtol::Float64=1e-5)
    num = 0.0
    den = 0.0
    (X <= 0 || abs(Npar) >= 1) && return num, den
    fnorm = μ / (4π * besselkx(2, μ))
    dn = 1 - Npar^2
    for n in 1:nmax
        A = n^2 * Y^2 + Npar^2 - 1
        A > 0 || continue
        upc = n * Y * Npar / dn
        apar = sqrt(A) / dn
        aperp = sqrt(A / dn)
        function integrand(φ)
            sφ, cφ = sincos(φ)
            upar = upc + apar * cφ
            uperp = aperp * sφ
            γ = n * Y + Npar * upar
            u2 = upar^2 + uperp^2
            dl = hypot(apar * sφ, aperp * cφ)
            grad = hypot(uperp / γ, upar / γ - Npar)
            wgt = 2π * uperp * uperp^2 * μ * fnorm * exp(-μ * (γ - 1)) / γ * dl / grad
            T = bessel_matrix(n, Nperp * uperp / Y, upar, uperp)
            W = wgt * real(dot(e, T * e))
            χu, χλ, λ = chi_derivatives(sf, upar, uperp, b)
            u = sqrt(u2)
            dχ = χu * γ / u + χλ * 2 * upar / (u2^2 * b) * (γ * upar - Npar * u2)
            return [W * dχ, W]
        end
        I, _ = quadgk(integrand, 0.0, π; rtol, atol=0.0, order=11)
        num += I[1]
        den += I[2]
    end
    return num, den
end

struct CurrentDriveTable
    ρ::Vector{Float64}
    sf::Vector{SpitzerFunction}
    G::Vector{Float64}      # F⟨1/R²⟩/⟨B⟩ [1/m] per surface, see `toroidal_factor`
    Rinv::Vector{Float64}   # ⟨1/R⟩ [1/m] per surface
end

"""
    CurrentDriveTable(m::PlasmaModel, Zeff; nρ=40, model=:full)

Responses on `nρ` surfaces. `model` is `:linliu` (the separable model of the Fortran,
`ncdroutine=1`), `:linliu_mc` (the same with the variational momentum-conserving Spitzer
function, `ncdroutine=2`), `:lorentz2d` (exact 2-D solution of the high-velocity operator,
`ncdroutine=3`) or `:full` (full linearized collision operator at the surface temperature,
`ncdroutine=4`)
"""
function CurrentDriveTable(m::PlasmaModel, Zeff::Real; nρ::Int=40, model::Symbol=:full)
    model in (:linliu, :linliu_mc, :lorentz2d, :full) || throw(ArgumentError("unknown current-drive model $model"))
    ρs = collect(range(0.04, 0.98; length=nρ))
    sfs = map(ρs) do ρ
        fs = FluxSurface(m, ρ)
        μ = 510.99895 / max(temperature(m, ρ), 1e-3)
        model == :linliu ? linliu_response(fs, Zeff, μ) :
        model == :linliu_mc ? linliu_response(fs, Zeff, μ; momentum_conservation=true) :
        model == :lorentz2d ? SpitzerFunction(fs, Zeff; nλ=300) :
        full_operator_response(fs, Zeff, μ)
    end
    F = m.R_axis * abs(B_cyl(m, m.R_axis, m.Z_axis)[2])      # R B_φ (flux function)
    G = [F * flux_average(sf.fs, 1 ./ sf.fs.R .^ 2) / flux_average(sf.fs, sf.fs.B) for sf in sfs]
    Rinv = [flux_average(sf.fs, 1 ./ sf.fs.R) for sf in sfs]
    return CurrentDriveTable(ρs, sfs, G, Rinv)
end

"""
    toroidal_factor(table::CurrentDriveTable, ρ)

`G = F ⟨1/R²⟩ / ⟨B⟩` [1/m] on the surface (`F = R B_φ`): the toroidal current is
`I = ∫ (⟨j∥⟩/⟨B⟩) dΨ_tor` with `dΨ_tor = ⟨B·∇φ⟩ dV / 2π = F ⟨1/R²⟩ dV / 2π`, so
`dI/dV = ⟨j∥⟩ G / 2π`, and the toroidal current density is `j_tor = ⟨j∥⟩ G / ⟨1/R⟩`
(Marushchenko et al. 2011, Appendix, Eqs. A1 and A9); `G ≈ 1/R₀` near the axis
"""
toroidal_factor(table::CurrentDriveTable, ρ::Real) = table_interp(table, table.G, ρ)

"""
    jtor_factor(table::CurrentDriveTable, ρ)

`G / ⟨1/R⟩`: converts the flux-surface-averaged `⟨j∥⟩` to the toroidal current density
"""
jtor_factor(table::CurrentDriveTable, ρ::Real) = table_interp(table, table.G ./ table.Rinv, ρ)

function table_interp(table::CurrentDriveTable, v::AbstractVector, ρ::Real)
    ρs = table.ρ
    k = clamp(searchsortedlast(ρs, ρ), 1, length(ρs) - 1)
    t = clamp((ρ - ρs[k]) / (ρs[k+1] - ρs[k]), 0.0, 1.0)
    return (1 - t) * v[k] + t * v[k+1]
end

"""
    cd_model(ncdroutine)

Current-drive response model for a `ncdroutine` value: 0/1 → `:linliu`, 2 → `:linliu_mc`
(the Fortran's models), 3 → `:lorentz2d`, 4 → `:full` (the exact solvers)
"""
cd_model(ncdroutine::Int) = ncdroutine <= 1 ? :linliu : ncdroutine == 2 ? :linliu_mc : ncdroutine == 3 ? :lorentz2d : :full

"""
    cd_efficiency(table::CurrentDriveTable, st, w, N; nmax=3)

Local `j∥/P_abs` at a plasma state, interpolated linearly in rho between the
two neighbouring tabulated surfaces
"""
function cd_efficiency(table::CurrentDriveTable, st, w::WaveParams, N::AbstractVector; nmax::Int=3, warm::Bool=true)
    ρ = sqrt(max(st.ψn, 0.0))
    ρ < 1 || return 0.0
    ρs = table.ρ
    if ρ <= ρs[1] || ρ >= ρs[end]
        k = ρ <= ρs[1] ? 1 : length(ρs)
        sf = table.sf[k]
        return cd_efficiency(sf, max(st.Bmag / sf.fs.Bmin, 1.0), st, w, N; nmax, warm)
    end
    k = searchsortedlast(ρs, ρ)
    t = (ρ - ρs[k]) / (ρs[k+1] - ρs[k])
    η1 = cd_efficiency(table.sf[k], max(st.Bmag / table.sf[k].fs.Bmin, 1.0), st, w, N; nmax, warm)
    η2 = cd_efficiency(table.sf[k+1], max(st.Bmag / table.sf[k+1].fs.Bmin, 1.0), st, w, N; nmax, warm)
    return (1 - t) * η1 + t * η2
end
