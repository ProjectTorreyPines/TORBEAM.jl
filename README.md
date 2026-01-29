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

The function `torbeam` that runs TORBEAM for all launchers with non-zero power for the current time point. Requires the environment variable $TORBEAM_DIR to be pointing to the `bin` folder with the TORBEAM executables.

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

## Usage instructions for Omega

Load the TORBEAM module with `module load torbeam` before using it.