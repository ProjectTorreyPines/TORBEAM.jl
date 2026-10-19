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

The function `run_torbeam` runs TORBEAM for all launchers with non-zero power for the current time point.

Outputs are stored in the `waves` and the `core_sources` IDS. 

## TorbeamParams

    Base.@kwdef struct TorbeamParams
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

## Fortran interface

`run_torbeam(dd, params)` calls the `beam` routine of `libtorbeamB.so`, found through
`TORBEAM_DIR` (`$TORBEAM_DIR/../lib/libtorbeamB.so`); `TORBEAM.fortran_available()`
tells whether the library can be found. The run is split into three steps:

- `equilibrium_inputs(dd)` assembles the launcher-independent part (`EquilibriumInputs`)
  and `beam_inputs(dd, ibeam, params, eq)` the per-launcher `BeamInputs` from IMAS,
- `run_beam(inputs, params)` runs one launcher and returns the raw `BeamOutputs`
  (the `rhoresult` scalars, the rays, the deposition profiles),
- `run_torbeam(dd, params)` drives all launchers and writes `waves` / `core_sources`.

Note on the `rhoresult` scalars (TORBEAM manual, 0-based numbering): `rhoresult(12)` is
the total driven current in kA, `rhoresult(13)` the absorbed power in MW and
`rhoresult(19)` the exit flag. Versions up to v1.0.5 read the current as the power
and a rho value as the current, so `waves.coherent_wave[].global_quantities.power`,
`current_tor` and the `core_sources` totals derived from them were wrong; they are
now `1e6 * rhoresult[14]` W and `1e3 * rhoresult[13]` A (1-based).

## Tests and golden data

`test/data/<case>.json` are trimmed `dd`s (equilibrium, core profiles, six EC beam
variants each for DIII-D 110 GHz X2/O2 and ITER 138 GHz O1) and `test/goldens/<case>.json`
the raw `BeamOutputs` the Fortran library produced for them, rounded to 10 significant
digits, with the wall time of each `beam` call.

`Pkg.test()` always checks the Julia-side input assembly against the goldens. When
`TORBEAM_DIR` points at the library it also checks `run_beam` (`rtol = 1e-6`) and the
end-to-end `run_torbeam`, and prints the measured wall time per beam next to the
recorded one (no test fails on speed). The Fortran library is only available on omega,
so these tests are not part of CI; run them there:

    module load torbeam/gcc11.x
    julia --project=. -e 'using Pkg; Pkg.test()'

Wall time per `beam` call recorded in the goldens (omega14, AMD EPYC 7502, libtorbeam
`c8/202505/build-generic`, 2026-10-19; a second run on the same node reproduced the
entries above 2 s within 7% and the totals within 2%):

| DIII-D 110 GHz             | [s]   | ITER 138 GHz               | [s]   |
|----------------------------|------:|----------------------------|------:|
| ec_pol20                   |  0.45 | ec_pol63                   |  26.5 |
| ec_pol20_focus_flipped     |  0.96 | ec_pol63_focus_flipped     |  29.9 |
| ec_pol10_tor15             |  2.91 | ec_pol60_tor20             |   9.4 |
| ec_pol30_tor-15            |  8.33 | ec_pol66_tor-20            | 106.6 |
| ec_pol30_tor15_half_power  |  4.12 | ec_pol63_tor20_half_power  |  12.5 |
| ec_pol30_tor15_Omode       |  2.39 | ec_pol66                   |  52.2 |
| total                      | 19.2  | total                      | 237.1 |

To regenerate the goldens on omega (reuses `test/data`, needs an environment with JSON
and this package `dev`ed):

    module load torbeam/gcc11.x
    julia --project=<env> test/goldens/generate.jl            # both cases, or: ... generate.jl D3D

With `--build` the `dd`s are rebuilt with FUSE first (the environment then also needs FUSE):

    julia --project=<env> test/goldens/generate.jl --build

## Usage instructions for Omega

Load the TORBEAM module with `module load torbeam` before using it.
