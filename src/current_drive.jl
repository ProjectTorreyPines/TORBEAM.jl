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
#     -ν_s u ∂χ/∂u + (ν_d/2) ⟨L_ξ⟩_b χ = -u ⟨ξ⟩_b ,
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
end

"""
    SpitzerFunction(fs::FluxSurface, Zeff; nu=300, nλ=200, umax=1.5, spitzer=nothing, μ=Inf)

Solve the bounce-averaged adjoint equation (see the file header) on the surface.

With `spitzer::SpitzerFunction1D` (the uniform-plasma Spitzer function with the
full linearized collision operator) and `μ = mc²/Te`, the solution is rescaled by
`spitzer_ratio(spitzer, u sqrt(μ/2))`, so that its u-dependence is that of the
full operator (electron-electron momentum conservation, energy diffusion, exact
thermal rates) while the pitch-angle/trapping structure is that of the
bounce-averaged Lorentz model — exact in a uniform field.
"""
function SpitzerFunction(fs::FluxSurface, Zeff::Real; nu::Int=300, nλ::Int=200, umax::Float64=1.5, spitzer=nothing, μ::Real=Inf)
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
        # implicit step: (χ_new - χ_old)/du = [ (νd/2) A χ_new + u ξbar ] / (νs u)
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
            rhs[k] = χ[iu-1, k] + c * u * ξbar[k]
        end
        χ[iu, 1:n] = tridiagonal_solve(lower, diag, upper, rhs)
    end
    if spitzer !== nothing
        @assert isfinite(μ) "the Spitzer rescaling needs μ = mc²/Te"
        for iu in 2:nu
            χ[iu, :] .*= spitzer_ratio(spitzer, us[iu] * sqrt(μ / 2))
        end
    end
    # spline on the (u, t) grid; evaluation maps λ -> t
    spl = Interpolations.cubic_spline_interpolation((us, ts), χ; extrapolation_bc=Interpolations.Line())
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
    cd_efficiency(sf::SpitzerFunction, b, st, w::WaveParams, N; nmax=3)

Local driven current per absorbed power, `j∥/P_abs` [A/W per m² → A m / W],
positive along B, at a point with `b = B/Bmin` on the surface of `sf`:
`-(e/(m_e c ν₀)) ⟨d·∇_u χ̂⟩_W` with ν₀ from the local density and the Coulomb logarithm
"""
function cd_efficiency(sf::SpitzerFunction, b::Real, st, w::WaveParams, N::AbstractVector; nmax::Int=3)
    X, Y = plasma_XY(st, w)
    (X <= 0 || st.Te <= 0) && return 0.0
    bhat = collect(st.B) ./ st.Bmag
    Npar = dot(N, bhat)
    Nperp = norm(N .- Npar .* bhat)
    μ = 510.99895 / st.Te
    e = cold_polarization(X, Y, Nperp, Npar)
    num, den = cd_resonance(X, Y, Nperp, Npar, μ, e, sf, b; nmax)
    den > 0 || return 0.0
    lnΛ = coulomb_log(st.ne, st.Te)
    ν0 = st.ne * e_charge^4 * lnΛ / (4π * ε_0^2 * m_e^2 * c_light^3)
    return -e_charge / (m_e * c_light * ν0) * num / den
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

"""
    CurrentDriveTable(m::PlasmaModel, Zeff; nρ=25)

Spitzer functions on `nρ` flux surfaces ρ_pol ∈ (0, 1), for interpolation along the ray
"""
struct CurrentDriveTable
    ρ::Vector{Float64}
    sf::Vector{SpitzerFunction}
end

"""
    CurrentDriveTable(m::PlasmaModel, Zeff; nρ=40, nλ=300, full_operator=true)

With `full_operator` the responses carry the u-dependence of the full linearized
collision operator (see `SpitzerFunction`), using the surface temperature
"""
function CurrentDriveTable(m::PlasmaModel, Zeff::Real; nρ::Int=40, nλ::Int=300, full_operator::Bool=true)
    ρs = collect(range(0.04, 0.98; length=nρ))
    spitzer = full_operator ? SpitzerFunction1D(Zeff) : nothing
    sfs = [SpitzerFunction(FluxSurface(m, ρ), Zeff; nλ, spitzer, μ=510.99895 / max(temperature(m, ρ), 1e-3)) for ρ in ρs]
    return CurrentDriveTable(ρs, sfs)
end

"""
    cd_efficiency(table::CurrentDriveTable, st, w, N; nmax=3)

Local `j∥/P_abs` at a plasma state, interpolated linearly in rho between the
two neighbouring tabulated surfaces
"""
function cd_efficiency(table::CurrentDriveTable, st, w::WaveParams, N::AbstractVector; nmax::Int=3)
    ρ = sqrt(max(st.ψn, 0.0))
    ρ < 1 || return 0.0
    ρs = table.ρ
    if ρ <= ρs[1] || ρ >= ρs[end]
        k = ρ <= ρs[1] ? 1 : length(ρs)
        sf = table.sf[k]
        return cd_efficiency(sf, max(st.Bmag / sf.fs.Bmin, 1.0), st, w, N; nmax)
    end
    k = searchsortedlast(ρs, ρ)
    t = (ρ - ρs[k]) / (ρs[k+1] - ρs[k])
    η1 = cd_efficiency(table.sf[k], max(st.Bmag / table.sf[k].fs.Bmin, 1.0), st, w, N; nmax)
    η2 = cd_efficiency(table.sf[k+1], max(st.Bmag / table.sf[k+1].fs.Bmin, 1.0), st, w, N; nmax)
    return (1 - t) * η1 + t * η2
end
