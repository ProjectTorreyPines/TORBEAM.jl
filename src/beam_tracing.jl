# Paraxial beam tracing (Pereverzev 1998; Poli, Peeters & Pereverzev 2001).
#
# The beam is a complex eikonal E ∝ exp(i k0 (s + iφ)): along the central ray
# x(τ), N(τ) the second derivatives of s and φ form the symmetric matrices S
# (wavefront curvature) and Φ (amplitude, i.e. width). With M = S + iΦ and the
# Hessian blocks of the dispersion function H, the central ray and M obey
#
#     dx/dτ = ∂H/∂N,  dN/dτ = -∂H/∂x,
#     dM/dτ = -(H_xx + H_xN M + M H_Nx + M H_NN M).
#
# Here everything is integrated in arclength s along the central ray
# (ds = |∂H/∂N| dτ). State vector: x (3), N (3), Re M and Im M (6 + 6, upper
# triangle). SI units; M in 1/m.

using OrdinaryDiffEqTsit5
import LinearAlgebra: dot, norm, cross, eigen, Symmetric, I

const NSTATE = 18
const M_IDX = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))

function pack_M!(u::AbstractVector, M::AbstractMatrix)
    for (k, (i, j)) in enumerate(M_IDX)
        u[6+k] = real(M[i, j])
        u[12+k] = imag(M[i, j])
    end
    return u
end

function unpack_M(u::AbstractVector)
    M = zeros(complex(eltype(u)), 3, 3)
    for (k, (i, j)) in enumerate(M_IDX)
        M[i, j] = M[j, i] = complex(u[6+k], u[12+k])
    end
    return M
end

"""
    Launch

Launch conditions of a beam: wave, position `x0` [m], unit direction `N0`,
local frame (`eh` horizontal, `ev` "vertical", both ⟂ to `N0`), 1/e widths
`wh`/`wv` [m] and wavefront curvature radii `Rh`/`Rv` [m] (TORBEAM sign:
positive = converging), launched power [W]
"""
struct Launch
    wave::WaveParams
    x0::Vector{Float64}
    N0::Vector{Float64}
    eh::Vector{Float64}
    ev::Vector{Float64}
    wh::Float64
    wv::Float64
    Rh::Float64
    Rv::Float64
    power::Float64
end

"""
    Launch(inputs::BeamInputs)

Launch conditions from TORBEAM's `floatinbeam`/`intinbeam` (cgs, degrees, MW):
`xpoldeg > 0` injects downwards, `xtordeg > 0` towards -y (clockwise from above)
"""
function Launch(inputs::BeamInputs)
    fi = inputs.floatinbeam
    wave = WaveParams(fi[1], inputs.intinbeam[3])
    θtor = deg2rad(fi[2])
    θpol = deg2rad(fi[3])
    x0 = [fi[4], fi[5], fi[6]] ./ 100
    N0 = [-cos(θpol) * cos(θtor), -cos(θpol) * sin(θtor), -sin(θpol)]
    eh = cross([0.0, 0.0, 1.0], N0)
    eh ./= norm(eh)
    ev = cross(N0, eh)
    return Launch(wave, x0, N0, eh, ev, fi[22] / 100, fi[23] / 100, fi[20] / 100, fi[21] / 100, fi[24] * 1e6)
end

"""
    initial_state(l::Launch)

State vector at the launch point: Gaussian beam with the given widths and
curvature radii. `Φ = 2/(k0 w²)` gives the 1/e amplitude radius `w`; `S = -1/R`
because a positive (converging) TORBEAM radius means the phase front is
centred ahead of the launch point.
"""
function initial_state(l::Launch)
    k0 = l.wave.k0
    Ph = l.eh * l.eh'
    Pv = l.ev * l.ev'
    S = -Ph / l.Rh - Pv / l.Rv
    Φ = 2 / (k0 * l.wh^2) * Ph + 2 / (k0 * l.wv^2) * Pv
    u0 = zeros(NSTATE)
    u0[1:3] = l.x0
    u0[4:6] = l.N0
    pack_M!(u0, S + im * Φ)
    return u0
end

struct BeamTracer{PM<:PlasmaModel}
    model::PM
    wave::WaveParams
    res::Any   # ForwardDiff DiffResults.HessianResult cache
end
BeamTracer(model::PlasmaModel, wave::WaveParams) = BeamTracer(model, wave, ForwardDiff.DiffResults.HessianResult(zeros(6)))

"""
    beam_rhs!(du, u, tracer::BeamTracer, s)

Right-hand side of the ray + Riccati equations in arclength
"""
function beam_rhs!(du, u, tracer::BeamTracer, s)
    res = dispersion_derivatives!(tracer.res, tracer.model, tracer.wave, @view u[1:6])
    g = ForwardDiff.DiffResults.gradient(res)
    Hh = ForwardDiff.DiffResults.hessian(res)
    gx = @view g[1:3]
    gN = @view g[4:6]
    vnorm = norm(gN)
    du[1:3] = gN / vnorm
    du[4:6] = -gx / vnorm
    M = unpack_M(u)
    Hxx = @view Hh[1:3, 1:3]
    HxN = @view Hh[1:3, 4:6]
    HNN = @view Hh[4:6, 4:6]
    dM = -(Hxx + HxN * M + M * HxN' + M * HNN * M) / vnorm
    pack_M!(du, dM)
    return du
end

"""
    BeamSolution

Result of `trace_beam`: the `OrdinaryDiffEq` solution (dense in arclength),
the launch, the total arclength and the exit reason (`:rhostop`, `:grid`, `:length`)
"""
struct BeamSolution{S}
    sol::S
    launch::Launch
    length::Float64
    exit::Symbol
end

"""
    trace_beam(m::PlasmaModel, l::Launch; rhostop=0.96, reltol=1e-7, abstol=1e-7, smax=20.0)

Integrate the central ray and the beam matrix from the launch point until the
beam, having entered the plasma, reaches `rhostop` on its way out, leaves the
equilibrium grid, or exceeds `smax` [m] of arclength.
"""
function trace_beam(m::PlasmaModel, l::Launch; rhostop::Float64=0.96, reltol::Float64=1e-7, abstol::Float64=1e-7, smax::Float64=20.0)
    tracer = BeamTracer(m, l.wave)
    u0 = initial_state(l)
    exit = Ref(:length)

    ρ_of(u) = rho_pol(m, hypot(u[1], u[2]), u[3])
    inside = Ref(false)
    # stop when rho crosses rhostop upwards after the beam has been inside
    function track_inside(u, s, integ)
        inside[] || (inside[] = ρ_of(u) < rhostop)
        return false
    end
    rho_cb = ContinuousCallback((u, s, integ) -> (inside[] ? ρ_of(u) - rhostop : -1.0),
        integ -> (exit[] = :rhostop; terminate!(integ)), nothing)
    inside_cb = DiscreteCallback(track_inside, integ -> nothing)
    grid_cb = DiscreteCallback((u, s, integ) -> begin
            R = hypot(u[1], u[2])
            !(m.R[1] <= R <= m.R[end] && m.Z[1] <= u[3] <= m.Z[end])
        end,
        integ -> (exit[] = :grid; terminate!(integ)))

    prob = ODEProblem(beam_rhs!, u0, (0.0, smax), tracer)
    sol = solve(prob, Tsit5(); reltol, abstol, callback=CallbackSet(inside_cb, rho_cb, grid_cb), dtmax=0.02)
    return BeamSolution(sol, l, sol.t[end], exit[])
end

"""
    beam_widths(b::BeamSolution, s; perp=:v)

Beam cross-section at arclength `s`: central ray position `x` [m] and the 1/e
amplitude half-widths [m] along three directions ⟂ to the group velocity
(`perp=:v`, default) or to the wave vector (`perp=:N`):

  - `wh` along `eh`, the horizontal direction ("left/right rays": the cut by the horizontal plane),
  - `wv` along `ev = v × eh`,
  - `wp` along `ep`, the direction in the poloidal (R-Z) plane ("upper/lower rays": the cut by the poloidal plane).
"""
function beam_widths(b::BeamSolution, s::Real; perp::Symbol=:v)
    u = b.sol(s)
    k0 = b.launch.wave.k0
    v = perp == :N ? u[4:6] : b.sol(s, Val{1})[1:3]
    v ./= norm(v)
    eh = cross([0.0, 0.0, 1.0], v)
    eh ./= norm(eh)
    ev = cross(v, eh)
    R = hypot(u[1], u[2])
    eR = [u[1] / R, u[2] / R, 0.0]
    ep = [0.0, 0.0, 1.0] * dot(v, eR) - eR * v[3]
    ep ./= norm(ep)
    Φ = imag(unpack_M(u))
    width(e) = sqrt(2 / (k0 * dot(e, Φ * e)))
    return (; x=u[1:3], wh=width(eh), wv=width(ev), wp=width(ep), eh, ev, ep)
end
