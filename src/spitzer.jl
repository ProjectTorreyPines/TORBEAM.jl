# Spitzer-Härm response function of a uniform, non-relativistic plasma with the
# full linearized electron collision operator (electron-electron test-particle
# and field-particle parts, electron-ion pitch-angle scattering with Z).
#
# The adjoint response for the parallel current, C†χ = -v ξ, is — because the
# linearized operator is self-adjoint under the f_M-weighted inner product —
# the forward Spitzer problem C(f_M χ) = -f_M v ξ: χ = ξ D(v) is the Spitzer
# function. In units v_T = sqrt(2T/m) = 1, ν̂ = n e⁴ lnΛ/(4π ε₀² m² v_T³) = 1,
# n = 1 (so f_M = π^{-3/2} e^{-x²}), with x = v/v_T:
#
#   test particle (Helander & Sigmar):
#     C_t(f₁) = -(ν_D^ee + ν_D^ei) f₁ + x^{-2} ∂_x[ x³ ( ½ν_s^ee f₁ + ½ν_∥^ee x ∂_x f₁ ) ],   f₁ = f_M D
#     ν_D^ee = (φ - G)/x³, ν_s^ee = 4G/x, ν_∥^ee = 2G/x³, ν_D^ei = Z/x³,
#     φ = erf(x), G = (φ - xφ')/(2x²) (Chandrasekhar).
#   field particle, from the l = 1 Rosenbluth potentials H₁ = ξ h, G₁ = ξ k of f₁:
#     C_f(f₁)/f_M = -x [ ψ'' + 4ψ'/x - 2(ψ + xψ') ],  ψ = h/x + x q',  q = k/x,
#     h(x) = 2π ∫ x'² f₁(x') I_h(x, x') dx',  k(x) = 2π ∫ x'² f₁(x') I_k(x, x') dx',
#     I_h = ∫₋₁¹ μ (a - bμ)^{-1/2} dμ,  I_k = ∫₋₁¹ μ (a - bμ)^{1/2} dμ,  a = x² + x'², b = 2xx'.
#
# The classic check is the conductivity ratio to the Lorentz gas,
# σ/σ_L = γ_E(Z) = 0.5816, 0.6833, 0.7849, 0.9252 for Z = 1, 2, 4, 16.

import SpecialFunctions: erf, gamma

chandrasekhar(x::Real) = (erf(x) - x * 2 / sqrt(π) * exp(-x^2)) / (2x^2)

"""
    SpitzerFunction1D(Z; n=600, xmax=7.0, sink=0.0)

Spitzer function `D(x)` (χ = ξ D, x = v/v_T) on a uniform grid, the Lorentz-gas
high-velocity solution `x⁴/(5+Z)` it tends to, and the conductivity ratio `γ_E`;
with `sink = f_tr/f_c` the l = 1 problem of a flux surface with trapped particles
"""
struct SpitzerFunction1D
    Z::Float64
    x::Vector{Float64}
    D::Vector{Float64}
    γE::Float64
    ee_residual::Float64   # relative residual of e-e momentum conservation on f_M v∥
end

"""
    derivative_matrices(x)

First and second derivative matrices on a uniform grid starting one step
after 0 (D(0) = 0 is the left neighbour of the first node), one-sided at the end
"""
function derivative_matrices(x::AbstractVector)
    n = length(x)
    dx = x[2] - x[1]
    D1 = zeros(n, n)
    D2 = zeros(n, n)
    for i in 1:n
        if i == 1
            D1[i, i+1] = 1 / (2dx)
            D2[i, i] = -2 / dx^2
            D2[i, i+1] = 1 / dx^2
        elseif i == n
            D1[i, i-1] = -1 / dx
            D1[i, i] = 1 / dx
            D2[i, i-2] = 1 / dx^2
            D2[i, i-1] = -2 / dx^2
            D2[i, i] = 1 / dx^2
        else
            D1[i, i-1] = -1 / (2dx)
            D1[i, i+1] = 1 / (2dx)
            D2[i, i-1] = 1 / dx^2
            D2[i, i] = -2 / dx^2
            D2[i, i+1] = 1 / dx^2
        end
    end
    return D1, D2
end

"""
    field_matrix(x)

Dense matrix of the l = 1 field-particle operator `C_f(f_M D ξ)/(f_M ξ)` acting on
`D` sampled on the uniform grid `x = v/v_T` (units ν̂, see the file header)
"""
function field_matrix(x::AbstractVector)
    n = length(x)
    dx = x[2] - x[1]
    fM = @. exp(-x^2) / π^1.5
    D1, D2 = derivative_matrices(x)
    H = zeros(n, n)
    K = zeros(n, n)
    for i in 1:n, j in 1:n
        a = x[i]^2 + x[j]^2
        b = 2 * x[i] * x[j]
        sp = x[i] + x[j]
        sm = abs(x[i] - x[j])
        Ih = (2a * (sp - sm) - (2 / 3) * (sp^3 - sm^3)) / b^2
        Ik = ((2a / 3) * (sp^3 - sm^3) - (2 / 5) * (sp^5 - sm^5)) / b^2
        H[i, j] = 2π * x[j]^2 * fM[j] * Ih * dx
        K[i, j] = 2π * x[j]^2 * fM[j] * Ik * dx
    end
    Xinv = Diagonal(1 ./ x)
    X = Diagonal(x)
    Q = Xinv * K                       # q = k/x
    Ψ = Xinv * H + X * (D1 * Q)        # ψ = h/x + x q'
    return -X * (D2 * Ψ + 4 * Xinv * (D1 * Ψ) - 2 * (Ψ + X * (D1 * Ψ)))
end

"""
    SpitzerFunction1D(Z; n=600, xmax=7.0)

Solve the l = 1 Spitzer-Härm problem with the full linearized operator
"""
function SpitzerFunction1D(Z::Real; n::Int=600, xmax::Float64=7.0, sink::Real=0.0)
    x = collect(range(0.0, xmax; length=n + 1))[2:end]     # x₁ = dx .. xmax; D(0) = 0
    dx = x[2] - x[1]
    fM = @. exp(-x^2) / π^1.5
    φ = erf.(x)
    G = chandrasekhar.(x)
    νDee = @. (φ - G) / x^3
    νs = @. 4G / x
    νpar = @. 2G / x^3
    νDei = @. Z / x^3
    D1, D2 = derivative_matrices(x)

    # test-particle operator on D (per f_M): -(νDee + νDei) D + (1/(x² f_M)) d/dx[ x³ (½νs f_M D + ½νpar x d(f_M D)/dx) ].
    # With (f_M D)' = f_M (D' - 2xD) the f_M factors cancel analytically:
    #   P = ½νs D + ½νpar x (D' - 2xD),  C_s = (x³P)'/x² - 2x² P
    P = Diagonal(0.5 .* νs) + Diagonal(@. 0.5 * νpar * x) * (D1 - Diagonal(2 .* x))
    T = -Diagonal(νDee .+ νDei) + Diagonal(@. 1 / x^2) * (D1 * (Diagonal(x .^ 3) * P)) - Diagonal(@. 2x^2) * P
    F = field_matrix(x)

    # e-e momentum conservation: the e-e part (T without the ion scattering, plus F) must
    # annihilate the shifted Maxwellian f₁ = f_M x ξ, i.e. (T + νDei + F) x = 0
    r = (T + Diagonal(νDei) + F) * x
    ee_residual = maximum(abs, r[4:n-3]) / maximum(abs, (νDee .* x)[4:n-3])   # interior nodes (boundary rows are one-sided)
    # optional trapped-particle momentum sink -(f_tr/f_c) ν_e D of the collisionless flux-surface
    # problem (Marushchenko et al. 2009 Eq. 3, Romé et al. 1998 Eq. A10), `sink = f_tr/f_c`
    A = T + F - sink * Diagonal(νDee .+ νDei)
    rhs = -x
    # at xmax impose the asymptotic power law D' = 4D/x (D -> x⁴/(5+Z) up to O(x⁻²))
    A[n, :] .= 0
    A[n, n-1] = -1 / dx
    A[n, n] = 1 / dx - 4 / x[n]
    rhs[n] = 0
    D = A \ rhs

    # conductivity ratio to the Lorentz gas (electron-ion scattering only)
    DL = @. x / νDei
    w = @. x^3 * fM
    γE = sum(w .* D) / sum(w .* DL)
    return SpitzerFunction1D(Float64(Z), x, D, γE, ee_residual)
end

"""
    spitzer_ratio_sink(sp::SpitzerFunction1D, x)

`D(x)` of the 1-D solution interpolated at `x`, with its asymptotic power law beyond the grid
"""
function spitzer_ratio_sink(sp::SpitzerFunction1D, x::Real)
    x >= sp.x[end] && return sp.D[end] * (x / sp.x[end])^4
    k = clamp(searchsortedlast(sp.x, x), 1, length(sp.x) - 1)
    t = (x - sp.x[k]) / (sp.x[k+1] - sp.x[k])
    return (1 - t) * sp.D[k] + t * sp.D[k+1]
end

"""
    spitzer_ratio(sp::SpitzerFunction1D, x)

`D(x) / (x⁴/(5+Z))`: the factor by which the full Spitzer function exceeds the
Lorentz-model high-velocity solution at `x = v/v_T` (-> 1 for x -> ∞)
"""
function spitzer_ratio(sp::SpitzerFunction1D, x::Real)
    x <= sp.x[1] && return spitzer_ratio(sp, sp.x[1])
    x >= sp.x[end] && return 1.0
    i = clamp(searchsortedlast(sp.x, x), 1, length(sp.x) - 1)
    t = (x - sp.x[i]) / (sp.x[i+1] - sp.x[i])
    D = sp.D[i] + t * (sp.D[i+1] - sp.D[i])
    return D / (x^4 / (5 + sp.Z))
end

# ---------------------------------------------------------------------------
# Variational Spitzer function with momentum conservation and a trapped-particle
# momentum sink (Romé et al., PPCF 40 (1998) 511, Appendix; weakly relativistic
# extension and explicit coefficients: Marushchenko, Beidler & Maassberg,
# Fusion Sci. Technol. 55 (2009) 180, Eqs. 5-8 and A2-A5). This is the model
# behind the Fortran's `ncdroutine=2`.
#
# The l = 1 Spitzer problem on a flux surface in the collisionless limit,
#     Ĉ₁^lin(K) - (f_tr/f_c) ν_e(u) K = ν_e0 (u/γ) F_eM ,      u = p/p_th,
# is solved with the trial function χ_a(u) = K/F_eM = (u/γ) Σ_{i=1}^4 d_i uⁱ by
# minimising the entropy-production functional with a Lagrange multiplier ζ for
# momentum conservation:
#     (M_ij + Ω_ij) d_j + (M_0i + Ω_0i) ζ = G_i   (i = 1..4),
#     (M_0j + Ω_0j) d_j                   = G_0 ,
# M_ij = M⁽⁰⁾ + M⁽¹⁾/μ (linearized e-e operator plus e-i pitch-angle scattering),
# Ω_ij = (f_tr/f_c)(ω⁽⁰⁾ + ω⁽¹⁾/μ) (the trapped-particle sink, ν_e = e-e + e-i
# pitch-angle rate), G_i = G⁽⁰⁾ + G⁽¹⁾/μ, μ = m c²/T_e. The μ⁰ parts are
# Hirshman's / Romé's non-relativistic coefficients. For f_c = 1 and μ → ∞ this
# is the classical Spitzer-Härm function in the normalisation of
# `SpitzerFunction1D` (u⁴/(5+Z) at high speed), with no free constant.
#
# The e-i (Z_eff) parts are the same operator moments in M and in Ω (M_ij,Z ≡ ω_ij,Z)
# and are evaluated analytically from F_eM ≈ π^{-3/2} e^{-u²} (1 + (u⁴/2 - 15/8)/μ) and
# the 1/γ of the test functions (one power for the i = 0 row, two otherwise):
#     ω⁽⁰⁾_ij,Z = Γ((n+2)/2),   ω⁽¹⁾_ij,Z = ½Γ((n+6)/2) - (15/8)Γ((n+2)/2) - c_i Γ((n+4)/2),
# n = i + j, c₀ = 1, c_{i≥1} = 2. The e-e part of row 0 is zero (momentum conservation,
# ∫ p∥ C^ee = 0); the other e-e parts are the published tables (Marushchenko 2009, A3/A4).
# The published M⁽¹⁾ rows 0 and 1 violate these identities and are not used.

"""
    variational_spitzer(fc, Zeff, μ)

Coefficients `d` (length 4) and the function `χ_a(u) = (u/γ) Σ dᵢ uⁱ`, `u = p/p_th`,
`γ = sqrt(1 + 2u²/μ)`, of the momentum-conserving Spitzer function on a surface with
circulating fraction `fc` (see the comment above)
"""
function variational_spitzer(fc::Real, Zeff::Real, μ::Real)
    sπ = sqrt(π)
    s2 = sqrt(2.0)
    r = (1 - fc) / fc                      # f_tr / f_c
    # e-i parts (analytic, see above), indices 0..4
    Z0(i, j) = Zeff * gamma((i + j + 2) / 2)
    Z1(i, j) = Zeff * (0.5gamma((i + j + 6) / 2) - 15 / 8 * gamma((i + j + 2) / 2) - (i == 0 ? 1 : 2) * gamma((i + j + 4) / 2))
    # e-e parts of M⁽⁰⁾, M⁽¹⁾ (A3), rows i ≥ 1 (row 0 vanishes by momentum conservation)
    M0 = zeros(5, 5)
    M1 = zeros(5, 5)
    M0[2, 2] = -104 / 15 + 151 / (15s2);          M1[2, 2] = 2071 / 105 - 197861 / (6720s2)
    M0[2, 3] = 4 / sπ - sπ;                       M1[2, 3] = -5 / sπ + sπ
    M0[2, 4] = -102 / 5 + 607 / (20s2);           M1[2, 4] = -1437 / 140 - 131497 / (8960s2)
    M0[2, 5] = 26 / sπ - 7sπ;                     M1[2, 5] = 52 / sπ - sπ / 4
    M0[3, 3] = s2;                                M1[3, 3] = -1105 / (256s2)
    M0[3, 4] = 6 / sπ;                            M1[3, 4] = -49 / (2sπ) + 69sπ / 8
    M0[3, 5] = 11 / s2;                           M1[3, 5] = -30245 / (1024s2)
    M0[4, 4] = -228 / 5 + 6147 / (80s2);          M1[4, 4] = 1161 / 14 - 550017 / (7168s2)
    M0[4, 5] = 45 / sπ - 9sπ / 4;                 M1[4, 5] = 79 / (4sπ) + 255sπ / 4
    M0[5, 5] = 157 / (2s2);                       M1[5, 5] = -2754977 / (4096s2)
    # e-e parts of ω⁽⁰⁾, ω⁽¹⁾ (A4); for i, j ≥ 1 they depend on i + j only
    W0 = zeros(5, 5)
    W1 = zeros(5, 5)
    W0[1, 2] = 1 / sπ;                            W1[1, 2] = -9 / (4sπ) + 3sπ / 4
    W0[1, 3] = 1 / s2;                            W1[1, 3] = 19 / (16s2)
    W0[1, 4] = 1 / sπ + sπ / 4;                   W1[1, 4] = 7 / (4sπ) + 7sπ / 8
    W0[1, 5] = 9 / (4s2);                         W1[1, 5] = 591 / (64s2)
    W0[2, 2] = W0[1, 3];                          W1[2, 2] = -17 / (16s2)
    W0[2, 3] = W0[1, 4];                          W1[2, 3] = -3 / (4sπ) + sπ / 8
    W0[2, 4] = W0[1, 5];                          W1[2, 4] = 131 / (64s2)
    W0[2, 5] = 5 / (2sπ) + 3sπ / 4;               W1[2, 5] = 41 / (8sπ) + 15sπ / 8
    W0[3, 3] = W0[2, 4];                          W1[3, 3] = W1[2, 4]
    W0[3, 4] = W0[2, 5];                          W1[3, 4] = W1[2, 5]
    W0[3, 5] = 115 / (16s2);                      W1[3, 5] = 7133 / (256s2)
    W0[4, 4] = W0[3, 5];                          W1[4, 4] = W1[3, 5]
    W0[4, 5] = 9 / sπ + 45sπ / 16;                W1[4, 5] = 203 / (4sπ) + 525sπ / 32
    W0[5, 5] = 1911 / (64s2);                     W1[5, 5] = 239853 / (1024s2)
    A = zeros(5, 5)
    for i in 1:5, j in i:5
        zz = Z0(i - 1, j - 1) + Z1(i - 1, j - 1) / μ
        A[i, j] = M0[i, j] + M1[i, j] / μ + zz + r * (W0[i, j] + W1[i, j] / μ + zz)
        A[j, i] = A[i, j]
    end
    # G (A5), i = 0..4
    G = [gamma((5 + i) / 2) + (-15 / 8 * gamma((5 + i) / 2) - 2gamma((7 + i) / 2) + gamma((9 + i) / 2) / 2) / μ for i in 0:4]
    # unknowns (d₁..d₄, ζ): rows 1..4 the i = 1..4 equations, row 5 the momentum constraint
    S = zeros(5, 5)
    rhs = zeros(5)
    for i in 1:4
        S[i, 1:4] = A[i+1, 2:5]
        S[i, 5] = A[1, i+1]
        rhs[i] = G[i+1]
    end
    S[5, 1:4] = A[1, 2:5]
    rhs[5] = G[1]
    d = (S \ rhs)[1:4]
    χ(u) = u / sqrt(1 + 2u^2 / μ) * (d[1] * u + d[2] * u^2 + d[3] * u^3 + d[4] * u^4)
    return d, χ
end
