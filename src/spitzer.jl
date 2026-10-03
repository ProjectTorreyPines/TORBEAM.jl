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

import SpecialFunctions: erf

chandrasekhar(x::Real) = (erf(x) - x * 2 / sqrt(π) * exp(-x^2)) / (2x^2)

"""
    SpitzerFunction1D

Spitzer function `D(x)` (χ = ξ D, x = v/v_T) on a uniform grid, the Lorentz-gas
high-velocity solution `x⁴/(5+Z)` it tends to, and the conductivity ratio `γ_E`
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
function SpitzerFunction1D(Z::Real; n::Int=600, xmax::Float64=7.0)
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
    A = T + F
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
