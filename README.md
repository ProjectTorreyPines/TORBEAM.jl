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
the published beam-tracing papers (no Fortran needed): cold-plasma paraxial
beam tracing (`src/dispersion.jl`, `src/beam_tracing.jl`), absorption from the
exactly relativistic anti-Hermitian dielectric tensor in the weak-damping
approximation, with the perpendicular index and polarization from the warm
(relativistic Hermitian) dispersion relation where the plasma is resonant
(`nabsroutine=1`; `nabsroutine=0` uses the cold polarization) (`src/absorption.jl`),
and deposition profiles from the beam's
Gaussian cross-section, with each part of the cross-section deposited where its
own path meets the resonance (`src/deposition.jl`). Against the Fortran it
reproduces the rays to a few mm and the deposition profiles (location, width,
shape) on both DIII-D-like (2 keV) and ITER-like (25 keV) cases. Current
drive (`src/current_drive.jl`) uses the adjoint method with a response function
solved numerically from the bounce-averaged adjoint Fokker-Planck equation on
each flux surface (relativistic test-particle collisions, Z_eff, trapping from
the real field variation). With `ncdroutine=1` this is the Lorentz-model
response (TORBEAM's Lin-Liu routine without momentum conservation, agreement
within ~20% for beams that drive significant current); with `ncdroutine=2`
(default) the response is that of the full linearized collision operator —
exact thermal rates with energy diffusion and the electron-electron
field-particle term (whose uniform-plasma limit, `src/spitzer.jl`,
reproduces the Spitzer-Härm conductivity ratios) — solved in the real trapped
geometry, where it reproduces the neoclassical conductivity (Sauter et al.
1999) within a few percent; it reproduces TORBEAM's momentum-conserving
current within ~16% on ITER-like and ~9% on DIII-D-like cases
(second-harmonic O-mode: right sign, ~40% high along with its weakly absorbed
power).

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

## Usage instructions for Omega

Load the TORBEAM module with `module load torbeam` before using it.