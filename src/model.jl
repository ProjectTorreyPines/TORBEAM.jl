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
    PlasmaModel(inputs::BeamInputs; edge_decay=0.02)

Equilibrium (psi, B on the R-Z grid) and kinetic profiles (ne, Te vs rho_pol)
interpolated with cubic splines, plus the magnetic axis found from the psi spline.

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
    ne::S1          # vs rho_pol [m^-3]
    Te::S1          # vs rho_pol [keV]
    ρ_edge::Float64
    ne_edge::Float64
    Te_edge::Float64
    edge_decay::Float64
    Zeff::Float64
    B0::Float64     # vacuum toroidal field at R0 [T]
    R0::Float64     # geometric axis [m]
    a::Float64      # minor radius [m]
end

function PlasmaModel(inputs::BeamInputs; edge_decay::Float64=0.02, nρ::Int=401)
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

    # profiles: prdata holds rho_pol, ne [1e19 m^-3], rho_pol, Te [keV]; the rho grid
    # is not uniform (uniform in psi), so resample on a uniform rho grid for the spline
    npsi = inputs.npsi
    prdata = inputs.prdata
    ρp = prdata[1:npsi]
    nep = prdata[npsi+1:2npsi] .* 1e19
    Tep = prdata[3npsi+1:4npsi]
    ρ = range(ρp[1], ρp[end]; length=nρ)
    ne = Interpolations.cubic_spline_interpolation(ρ, linear_resample(ρp, nep, ρ); extrapolation_bc=Interpolations.Line())
    Te = Interpolations.cubic_spline_interpolation(ρ, linear_resample(ρp, Tep, ρ); extrapolation_bc=Interpolations.Line())

    fi = inputs.floatinbeam
    return PlasmaModel(R, Z, ψ, BR, Bφ, BZ, ψ_axis, ψ_boundary, R_axis, Z_axis,
        ne, Te, ρp[end], nep[end], Tep[end], edge_decay, fi[35], fi[27], fi[25] / 100, fi[26] / 100)
end

"""
    linear_resample(x, y, xnew)

Piecewise-linear resampling of `y(x)` (x increasing) on the points `xnew`
"""
function linear_resample(x::AbstractVector, y::AbstractVector, xnew::AbstractVector)
    out = similar(y, length(xnew))
    for (k, xk) in enumerate(xnew)
        i = clamp(searchsortedlast(x, xk), 1, length(x) - 1)
        t = (xk - x[i]) / (x[i+1] - x[i])
        out[k] = y[i] + t * (y[i+1] - y[i])
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
    density(m::PlasmaModel, ρ)

Electron density [m^-3] vs rho_pol, exponentially decaying beyond the last profile point
"""
function density(m::PlasmaModel, ρ::Real)
    if ρ <= m.ρ_edge
        return m.ne(ρ)
    else
        return m.ne_edge * exp(-(ρ - m.ρ_edge) / m.edge_decay)
    end
end

"""
    temperature(m::PlasmaModel, ρ)

Electron temperature [keV] vs rho_pol, held at its edge value beyond the last profile point
"""
temperature(m::PlasmaModel, ρ::Real) = ρ <= m.ρ_edge ? m.Te(ρ) : m.Te_edge + zero(ρ)

"""
    state(m::PlasmaModel, x, y, z)

Local plasma state at the Cartesian point `(x, y, z)`: `B` [T], `Bmag` [T], `ne` [m^-3], `Te` [keV], `ρ`, `R`, `Z`
"""
function state(m::PlasmaModel, x::Real, y::Real, z::Real)
    R = hypot(x, y)
    ρ = rho_pol(m, R, z)
    B = B_cart(m, x, y, z)
    Bmag = sqrt(B[1]^2 + B[2]^2 + B[3]^2)
    return (; B, Bmag, ne=density(m, ρ), Te=temperature(m, ρ), ρ, R, Z=z)
end
