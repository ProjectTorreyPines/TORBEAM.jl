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

## Usage instructions for Omega

Load the TORBEAM module with `module load torbeam` before using it.