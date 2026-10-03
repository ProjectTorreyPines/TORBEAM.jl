# Plasma / equilibrium model for the Julia backend.
#
# Built from the same `BeamInputs` the Fortran library consumes (B components and
# psi on the rectangular R-Z grid, ne/Te vs rho_pol), so both backends see
# identical inputs. Everything is SI (m, T, m^-3) except Te in keV. All
# evaluation functions are generic in the coordinate type so ForwardDiff can
# differentiate through them.

import Interpolations
import ForwardDiff

const Spline2D = Interpolations.Extrapolation
const Spline1D = Interpolations.Extrapolation

"""
    PlasmaModel(inputs::BeamInputs; edge_decay=0.02, profile_interp=:cubic)

Equilibrium (psi, B on the R-Z grid) and kinetic profiles (ne, Te) interpolated
with cubic splines, plus the magnetic axis found from the psi spline.

The profiles are splined against the normalized flux `ψn = ρ_pol²` (uniform
grid), not ρ_pol: no square root is taken along the beam, so the derivatives
ForwardDiff propagates stay finite through the magnetic axis.

Beyond the last profile point `ρ_edge` (normally 1, the separatrix) the density
decays as `exp(-(ρ-ρ_edge)/edge_decay)` and Te is held at its edge value, so the
beam sees a smooth vacuum-plasma transition.
"""
struct PlasmaModel{RR<:AbstractRange,S2,S1}
    R::RR
    Z::RR
    ψ::S2
    BR::S2
    Bφ::S2
    BZ::S2
    ψ_axis::Float64
    ψ_boundary::Float64
    R_axis::Float64
    Z_axis::Float64
    ne::S1          # vs ψn = rho_pol² [m^-3]
    Te::S1          # vs ψn = rho_pol² [keV]
    ρ_edge::Float64
    ne_edge::Float64
    Te_edge::Float64
    edge_decay::Float64
    Zeff::Float64
    B0::Float64     # vacuum toroidal field at R0 [T]
    R0::Float64     # geometric axis [m]
    a::Float64      # minor radius [m]
end

function PlasmaModel(inputs::BeamInputs; edge_decay::Float64=0.02, nρ::Int=401, profile_interp::Symbol=:cubic)
    ni, nj = inputs.ni, inputs.nj
    eqdata = inputs.eqdata
    ψ_boundary = eqdata[1]
    Rv = eqdata[2:ni+1]
    Zv = eqdata[ni+2:ni+nj+1]
    block(b) = reshape(eqdata[1+ni+nj+b*ni*nj+1:1+ni+nj+(b+1)*ni*nj], ni, nj)
    BRg, Bφg, BZg, ψg = block(0), block(1), block(2), block(3)

    R = range(Rv[1], Rv[end]; length=ni)
    Z = range(Zv[1], Zv[end]; length=nj)
    @assert maximum(abs.(R .- Rv)) < 1e-9 * (Rv[end] - Rv[1]) "R grid must be uniform"
    @assert maximum(abs.(Z .- Zv)) < 1e-9 * (Zv[end] - Zv[1]) "Z grid must be uniform"
    spline2d(M) = Interpolations.cubic_spline_interpolation((R, Z), M; extrapolation_bc=Interpolations.Line())
    ψ, BR, Bφ, BZ = spline2d(ψg), spline2d(BRg), spline2d(Bφg), spline2d(BZg)

    # magnetic axis: psi is a minimum there when sgnm > 0 (psi_boundary > psi_axis)
    sgnm = inputs.floatinbeam[34]
    R_axis, Z_axis, ψ_axis = find_axis(ψ, R, Z, ψg, sgnm)

    # profiles: prdata holds rho_pol, ne [1e19 m^-3], rho_pol, Te [keV] on a
    # non-uniform rho grid; resample them (natural cubic spline, so the curvature
    # stays smooth) on a uniform grid in ψn = rho_pol² for the evaluation splines
    npsi = inputs.npsi
    prdata = inputs.prdata
    ρp = prdata[1:npsi]
    nep = prdata[npsi+1:2npsi] .* 1e19
    Tep = prdata[3npsi+1:4npsi]
    if profile_interp == :cubic
        ψn = range(ρp[1]^2, ρp[end]^2; length=nρ)
        ne = Interpolations.cubic_spline_interpolation(ψn, cubic_resample(ρp, nep, sqrt.(ψn)); extrapolation_bc=Interpolations.Line())
        Te = Interpolations.cubic_spline_interpolation(ψn, cubic_resample(ρp, Tep, sqrt.(ψn)); extrapolation_bc=Interpolations.Line())
    elseif profile_interp == :linear
        # piecewise linear in rho_pol on the input nodes (zero curvature inside segments)
        ne = Interpolations.extrapolate(Interpolations.interpolate((ρp .^ 2,), nep, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Line())
        Te = Interpolations.extrapolate(Interpolations.interpolate((ρp .^ 2,), Tep, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Line())
    else
        error("profile_interp must be :cubic or :linear")
    end

    fi = inputs.floatinbeam
    return PlasmaModel(R, Z, ψ, BR, Bφ, BZ, ψ_axis, ψ_boundary, R_axis, Z_axis,
        ne, Te, ρp[end], nep[end], Tep[end], edge_decay, fi[35], fi[27], fi[25] / 100, fi[26] / 100)
end

"""
    cubic_resample(x, y, xnew)

Natural cubic spline through `y(x)` (x increasing) evaluated at `xnew`
"""
function cubic_resample(x::AbstractVector, y::AbstractVector, xnew::AbstractVector)
    n = length(x)
    h = diff(x)
    # second derivatives M from the tridiagonal system (natural: M[1] = M[n] = 0)
    A = zeros(n, n)
    rhs = zeros(n)
    A[1, 1] = 1
    A[n, n] = 1
    for i in 2:n-1
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        rhs[i] = 6 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    end
    M = A \ rhs
    out = similar(y, length(xnew))
    for (k, xk) in enumerate(xnew)
        i = clamp(searchsortedlast(x, xk), 1, n - 1)
        t = x[i+1] - xk
        u = xk - x[i]
        out[k] = (M[i] * t^3 + M[i+1] * u^3) / (6h[i]) + (y[i] / h[i] - M[i] * h[i] / 6) * t + (y[i+1] / h[i] - M[i+1] * h[i] / 6) * u
    end
    return out
end

"""
    find_axis(ψ, R, Z, ψg, sgnm)

Magnetic axis as the extremum of the psi spline (minimum for `sgnm > 0`):
Newton iterations on ∇ψ = 0 from the extremal grid point.
"""
function find_axis(ψ, R::AbstractRange, Z::AbstractRange, ψg::AbstractMatrix, sgnm::Real)
    idx = sgnm > 0 ? argmin(ψg) : argmax(ψg)
    r, z = R[idx[1]], Z[idx[2]]
    for _ in 1:50
        g = Interpolations.gradient(ψ, r, z)
        H = Interpolations.hessian(ψ, r, z)
        δ = H \ g
        r -= δ[1]
        z -= δ[2]
        abs(δ[1]) + abs(δ[2]) < 1e-12 && break
    end
    @assert R[1] < r < R[end] && Z[1] < z < Z[end] "magnetic axis search left the grid: (R, Z) = ($r, $z)"
    return r, z, ψ(r, z)
end

# ---------------------------------------------------------------------------
# evaluation (generic in the coordinate type for ForwardDiff)
# ---------------------------------------------------------------------------

psi(m::PlasmaModel, R::Real, Z::Real) = m.ψ(R, Z)

"""
    rho_pol(m::PlasmaModel, R, Z)

Normalized poloidal flux coordinate sqrt((ψ-ψ_axis)/(ψ_boundary-ψ_axis)), clipped at 0
"""
function rho_pol(m::PlasmaModel, R::Real, Z::Real)
    ψn = (psi(m, R, Z) - m.ψ_axis) / (m.ψ_boundary - m.ψ_axis)
    return sqrt(max(ψn, zero(ψn)))
end

"""
    B_cyl(m::PlasmaModel, R, Z)

Magnetic field components `(B_R, B_φ, B_Z)` [T]
"""
B_cyl(m::PlasmaModel, R::Real, Z::Real) = (m.BR(R, Z), m.Bφ(R, Z), m.BZ(R, Z))

"""
    B_cart(m::PlasmaModel, x, y, z)

Magnetic field components `(B_x, B_y, B_z)` [T] at the Cartesian point `(x, y, z)`
"""
function B_cart(m::PlasmaModel, x::Real, y::Real, z::Real)
    R = hypot(x, y)
    cφ, sφ = x / R, y / R
    BR, Bφ, BZ = B_cyl(m, R, z)
    return (BR * cφ - Bφ * sφ, BR * sφ + Bφ * cφ, BZ)
end

"""
    psi_norm(m::PlasmaModel, R, Z)

Normalized poloidal flux (ψ-ψ_axis)/(ψ_boundary-ψ_axis), i.e. rho_pol²
"""
psi_norm(m::PlasmaModel, R::Real, Z::Real) = (psi(m, R, Z) - m.ψ_axis) / (m.ψ_boundary - m.ψ_axis)

"""
    density_ψn(m::PlasmaModel, ψn)

Electron density [m^-3] vs normalized flux, exponentially decaying (in rho_pol) beyond the last profile point
"""
function density_ψn(m::PlasmaModel, ψn::Real)
    if ψn <= m.ρ_edge^2
        return m.ne(ψn)
    else
        return m.ne_edge * exp(-(sqrt(ψn) - m.ρ_edge) / m.edge_decay)
    end
end

"""
    temperature_ψn(m::PlasmaModel, ψn)

Electron temperature [keV] vs normalized flux, held at its edge value beyond the last profile point
"""
temperature_ψn(m::PlasmaModel, ψn::Real) = ψn <= m.ρ_edge^2 ? m.Te(ψn) : m.Te_edge + zero(ψn)

"""
    density(m::PlasmaModel, ρ)

Electron density [m^-3] vs rho_pol
"""
density(m::PlasmaModel, ρ::Real) = density_ψn(m, ρ^2)

"""
    temperature(m::PlasmaModel, ρ)

Electron temperature [keV] vs rho_pol
"""
temperature(m::PlasmaModel, ρ::Real) = temperature_ψn(m, ρ^2)

"""
    state(m::PlasmaModel, x, y, z)

Local plasma state at the Cartesian point `(x, y, z)`: `B` [T], `Bmag` [T], `ne` [m^-3], `Te` [keV], `ψn` (= rho_pol²), `R`, `Z`
"""
function state(m::PlasmaModel, x::Real, y::Real, z::Real)
    R = hypot(x, y)
    ψn = psi_norm(m, R, z)
    B = B_cart(m, x, y, z)
    Bmag = sqrt(B[1]^2 + B[2]^2 + B[3]^2)
    return (; B, Bmag, ne=density_ψn(m, ψn), Te=temperature_ψn(m, ψn), ψn, R, Z=z)
end
