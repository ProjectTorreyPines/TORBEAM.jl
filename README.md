# TORBEAM.jl

Run the TORBEAM beam tracing code from Julia

    @article{poli2018torbeam,
    title={TORBEAM 2.0, a paraxial beam tracing code for electron-cyclotron beams in fusion plasmas for extended physics applications},
    author={Poli, Emanuele and Bock, A and Lochbrunner, M and Maj, Omar and Reich, M and Snicker, A and Stegmeir, Andreas and Volpe, F and Bertelli, Nicola and Bilato, Roberto and others},
    journal={Computer Physics Communications},
    volume={225},
    pages={36--46},
    year={2018},
    publisher={Elsevier}
    }

This package calls the FORTRAN API.

The function `torbeam` that runs TORBEAM for all launchers with non-zero power for the current time point.

Outputs are stored in the `waves` and the `core_sources` IDS. 

## TorbeamParams

    Base.@kwdef mutable struct TorbeamParams
        # switches
        npow::Int = 1             # Power absorption switch (1 = on, 0 = off)
        ncd::Int = 1              # Current drive calculation switch (1 = on, 0 = off)
        ncdroutine::Int = 2       # Current drive routine selection (0 = Curba, 1 = Lin-Liu, 2 = Lin-Liu + momentum conservation)
        nprofv::Int = 50          # Number of radial points for volume profile calculation
        noout::Int = 0            # Screen output switch (0 = output enabled, 1 = output disabled)
        nrela::Int = 1            # Relativity consideration in absorption (0 = weakly, 1 = fully relativistic)
        nmaxh::Int = 3            # Number of harmonics to consider (1 to 5)
        nabsroutine::Int = 1      # Absorption routine selection (0 = Westerhof, 1 = Farina)
        nastra::Int = 0           # Definition of driven current density (0 = Lin-Liu, 1 = ASTRA, 2 = JINTRAC)
        nprofcalc::Int = 1        # Deposition profile calculation method (0 = standard, 1 = Maj method)
        ncdharm::Int = 1          # Harmonic consideration in current drive efficiency (0 = lowest harmonic only, 1 = includes next harmonic)
        nrel::Int = 0             # Relativistic mass correction for reflectometry (1 = enabled, 0 = disabled)
        n_ray::Int = 5            # Number of rays used in beam tracing
        verbose::Bool = false     # Verbose mode for output debugging (true = enabled, false = disabled)

        # Float parameters
        xrtol::Float64 = 1e-07    # Required relative error tolerance
        xatol::Float64 = 1e-07    # Required absolute error tolerance
        xstep::Float64 = 2.0      # Integration step in vacuum (cm)
        rhostop::Float64 = 0.96   # Maximum value of the flux coordinate (rho) before stopping
        xzsrch::Float64 = 0.0     # Vertical position for searching the magnetic axis (default 0 cm)
    end

## Backends

`TorbeamParams(backend=:fortran)` (default) calls the `beam` routine of
`libtorbeamB.so`, found through `TORBEAM_DIR` (`$TORBEAM_DIR/../lib/libtorbeamB.so`).
`TORBEAM.fortran_available()` tells whether the library can be found.

`TorbeamParams(backend=:julia)` runs a pure-Julia implementation written from
the published papers (no Fortran needed; see [References](#references)):
cold-plasma paraxial beam tracing (`src/dispersion.jl`, `src/beam_tracing.jl`),
absorption from the exactly relativistic anti-Hermitian dielectric tensor in
the weak-damping approximation, with the perpendicular index and polarization
from the warm (relativistic Hermitian) dispersion relation where the plasma is
resonant (`src/absorption.jl`), deposition profiles from the beam's Gaussian
cross-section (`src/deposition.jl`) and current drive by the adjoint method
(`src/current_drive.jl`, `src/spitzer.jl`). Against the Fortran it reproduces
the rays to a few mm and the deposition profiles (location, width, shape) on
both DIII-D-like (2 keV) and ITER-like (25 keV) cases.

The switches keep their Fortran meaning, so the two backends can be compared
setting by setting, and the Julia backend adds higher-fidelity options on top:

| switch | value | model | agreement with the Fortran |
|---|---|---|---|
| `ncdroutine` | 1 | Lin-Liu et al. (2003): separable response χ = sgn(u∥) F(u) H(λ), slowing-down kept to its l = 1 moment (circulating fraction f_c), relativistic high-speed limit | ITER within 2 % (one beam 5 %), DIII-D X2 +9–13 %, near-perpendicular launches (near-cancelling currents) ×1.7–2.2 |
| | 2 (default) | the same with momentum conservation: the variational Spitzer function of Romé et al. (1998) with the trapped-particle momentum sink, as the non-relativistic enhancement over the high-speed limit (`variational_spitzer`) | ITER within 2.5 % (one beam 5 %), DIII-D X2 +8–12 %; the enhancement itself matches to 1 % on both |
| | 3 | exact 2-D (u, λ) solution of the bounce-averaged adjoint equation with the same relativistic high-velocity operator | validated against the Lorentz-gas conductivity 1 − f_t (exact) and the separable model in the uniform limit |
| | 4 | full linearized collision operator (exact thermal rates, energy diffusion, e–e field term, relativistic detailed balance) in the real trapped geometry | reproduces the Spitzer–Härm conductivity ratios and the neoclassical conductivity of Sauter et al. (1999) within a few %; 5–20 % below the Fortran's momentum-conserving currents |
| `nabsroutine` | 1 (default) | warm (relativistic Hermitian) N⊥ and polarization in the resonant layer, absorption from the complex root of the full relativistic dispersion relation, α = 2k₀ Im N⊥ (x̂·v̂) (Farina's WARMDISP route) | absorbed powers within 1 % (DIII-D O2 1.27 vs 1.25 MW), deposition medians to ≤ 0.005 on both machines |
| | 2 | the same N⊥ and polarization with the weak-damping absorption α = 2k₀κ, κ = −(e*ε^a e)/(v̂·∂λ/∂N) | 2–4× cheaper; identical on ITER, DIII-D O2 absorbs 12 % too much |
| | 0 | cold N⊥ and polarization | fast path |
| `nprofcalc` | 1 (default) | Poli et al. 2018 Eq. 14: each absorption step spread over the beam amplitude on the vertical plane through the step (the resonance taken as vertical), every sample on its own flux surface | ITER widths (16–84 %) within 10 %, peaks within 8 %; DIII-D medians to 0.001, widths 20–30 % wider |
| | 2 | the same on the local iso-Y surface (the actual resonance surface), shift limited to 2ξ for grazing crossings | narrower by 15–25 % than the Fortran |

The driven current is reported as the toroidal current density
j_tor = ⟨j∥⟩ F⟨1/R²⟩/(⟨B⟩⟨1/R⟩) and the total as ∫ (⟨j∥⟩/⟨B⟩) dΨ_tor, which is
how the Fortran's totals and profiles relate. Known differences left on
DIII-D: profile widths 20–30 % wider than the Fortran's (beam widths agree
to 5 %), X2 currents +10–13 %, the near-cancelling near-perpendicular currents
×2–3, second-harmonic O-mode current +34 %.

The run is split into three steps that a backend plugs into:

- `equilibrium_inputs(dd)` / `beam_inputs(dd, ibeam, params, eq)` assemble the
  TORBEAM input vectors (`BeamInputs`) from IMAS,
- `run_beam(inputs, params)` runs one launcher and returns the raw `BeamOutputs`,
- `run_torbeam(dd, params)` drives all launchers and writes `waves` / `core_sources`.

## Tests and golden data

`test/data/<case>.json` are trimmed `dd`s (equilibrium, core profiles, several EC
beam variants) and `test/goldens/<case>.json` the raw `BeamOutputs` the Fortran
library produced for them. `Pkg.test()` always checks the Julia-side input
assembly against the goldens, and additionally the Fortran backend when
`TORBEAM_DIR` points at the library (on omega: `module load torbeam`).

To regenerate the goldens (needs FUSE, run on omega):

    module load torbeam/gcc11.x
    julia --project=<env with FUSE, JSON and this package dev'ed> test/goldens/generate.jl

## References

The Julia backend is a clean-room implementation: it was written from the
published descriptions below and from first principles, and the Fortran source
was not consulted (the library is used only as a black box for the golden
comparisons). Where the backend offers the same reduced model as the Fortran
it follows the paper the Fortran cites; the higher-fidelity options are
validated against the classical results listed last.

Beam tracing and the TORBEAM models

- G. V. Pereverzev, *Beam tracing in inhomogeneous anisotropic plasmas*,
  Phys. Plasmas 5 (1998) 3529 — paraxial WKB (complex-eikonal) beam tracing.
- E. Poli, A. G. Peeters, G. V. Pereverzev, *TORBEAM, a beam tracing code for
  electron-cyclotron waves in tokamak plasmas*, Comput. Phys. Commun. 136
  (2001) 90 — the 19 beam-tracing ODEs and the absorption equation.
- E. Poli et al., *TORBEAM 2.0, a paraxial beam tracing code for
  electron-cyclotron beams in fusion plasmas for extended physics
  applications*, Comput. Phys. Commun. 225 (2018) 36 — which absorption,
  current-drive and deposition models the switches select (Sections 4 and 5,
  Eq. 14 for the `nprofcalc=1` profile).

Dielectric tensor and absorption

- T. H. Stix, *Waves in Plasmas* (AIP, 1992) — cold tensor, harmonic (Bessel)
  matrix, resonance geometry; the relativistic Maxwellian tensors are derived
  from the gyro-orbit integrals in `src/absorption.jl`.

Current drive

- T. M. Antonsen and K. R. Chu, Phys. Fluids 25 (1982) 1295 — the adjoint
  (Green's-function) formulation of rf current drive.
- Y. R. Lin-Liu, V. S. Chan, R. Prater, *Electron cyclotron current drive
  efficiency in general tokamak geometry*, Phys. Plasmas 10 (2003) 4064
  (GA report A24257) — the separable response χ = sgn(u∥) F(u) H(λ) with the
  circulating fraction f_c (`ncdroutine=1`), the efficiency Eqs. 38–40 and
  the ⟨j∥B⟩/⟨B²⟩ current definition.
- M. Romé, V. Erckmann, U. Gasparino, N. Karulin, Plasma Phys. Control.
  Fusion 40 (1998) 511, Appendix — the variational (fifth-degree polynomial)
  Spitzer function with momentum conservation and the trapped-particle
  momentum sink (f_tr/f_c) ν_e.
- N. B. Marushchenko, C. D. Beidler, H. Maassberg, *Current drive
  calculations with an advanced adjoint approach*, Fusion Sci. Technol. 55
  (2009) 180 — the weakly relativistic (μ⁻¹) extension of that variational
  Spitzer function, with the explicit matrix coefficients (`ncdroutine=2`).
- N. B. Marushchenko et al., *Electron cyclotron current drive in low
  collisionality limit: on parallel momentum conservation*, Phys. Plasmas 18
  (2011) 032501 — the relation between the high-speed-limit, the
  momentum-conserving and the exact bounce-averaged solutions, and the
  toroidal-current conventions (Appendix).
- L. Spitzer and R. Härm, Phys. Rev. 89 (1953) 977 — conductivity ratios
  γ_E(Z) used to validate the uniform-field solver.
- O. Sauter, C. Angioni, Y. R. Lin-Liu, Phys. Plasmas 6 (1999) 2834 — the
  collisionless neoclassical conductivity used to validate the full-operator
  solver in the real trapped geometry.

## Usage instructions for Omega

Load the TORBEAM module with `module load torbeam` before using it.
