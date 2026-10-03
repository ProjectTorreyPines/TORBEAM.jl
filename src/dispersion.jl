# Cold-plasma dispersion function for the beam tracer.
#
# H(x, N) = N² - n²(x, θ) with n² the Appleton-Hartree refractive index of the
# selected mode (O: +, X: -) for the angle θ between N and B. H = 0 on the
# dispersion surface; the rays are its Hamiltonian flow. All derivatives come
# from ForwardDiff through `PlasmaModel`.

const e_charge = 1.602176634e-19
const m_e = 9.1093837015e-31
const ε_0 = 8.8541878128e-12
const c_light = 2.99792458e8

"""
    WaveParams(frequency, mode)

Wave frequency [Hz] and mode (`1` = O, `-1` = X, as `ec_launchers.beam[].mode`)
"""
struct WaveParams
    ω::Float64
    k0::Float64      # vacuum wavenumber ω/c [1/m]
    mode::Int
end
WaveParams(frequency::Real, mode::Integer) = WaveParams(2π * frequency, 2π * frequency / c_light, Int(mode))

"""
    plasma_XY(s, w::WaveParams)

Stix parameters `X = ω_pe²/ω²` and `Y = ω_ce/ω` from a plasma `state`
"""
function plasma_XY(s, w::WaveParams)
    X = s.ne * e_charge^2 / (ε_0 * m_e * w.ω^2)
    Y = e_charge * s.Bmag / (m_e * w.ω)
    return X, Y
end

"""
    refractive_index2(X, Y, cos2θ, mode)

Appleton-Hartree cold-plasma `n²` for the ordinary (`mode = 1`) or extraordinary (`-1`) wave
"""
function refractive_index2(X::Real, Y::Real, cos2θ::Real, mode::Integer)
    sin2θ = 1 - cos2θ
    disc = sqrt(0.25 * Y^4 * sin2θ^2 + (1 - X)^2 * Y^2 * cos2θ)
    den = 1 - X - 0.5 * Y^2 * sin2θ + mode * disc
    return 1 - X * (1 - X) / den
end

"""
    dispersion(m::PlasmaModel, w::WaveParams, u)

Dispersion function `H = N² - n²` at `u = [x, y, z, Nx, Ny, Nz]`
"""
function dispersion(m::PlasmaModel, w::WaveParams, u::AbstractVector)
    s = state(m, u[1], u[2], u[3])
    X, Y = plasma_XY(s, w)
    N2 = u[4]^2 + u[5]^2 + u[6]^2
    NB = u[4] * s.B[1] + u[5] * s.B[2] + u[6] * s.B[3]
    cos2θ = NB^2 / (N2 * s.Bmag^2)
    return N2 - refractive_index2(X, Y, cos2θ, w.mode)
end

"""
    dispersion_derivatives!(res, m, w, u)

Value, gradient (6) and Hessian (6x6) of the dispersion function at `u`,
stored in the `ForwardDiff.DiffResults.HessianResult` `res`
"""
function dispersion_derivatives!(res, m::PlasmaModel, w::WaveParams, u::AbstractVector)
    return ForwardDiff.hessian!(res, v -> dispersion(m, w, v), u)
end
