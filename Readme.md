# TORBEAM.jl instructions
The `TORBEAM,jl` module provides the function `torbeam` that runs TORBEAM for all launchers with non-zero power for the current time point. Outputs are stored in the `waves` and the `core_sources` IDS. 

# Usage instructions for Omega
It's important to load the TORBEAM module on omega before the code is run, e.g.:
```
module load torbeam
julia --project=<path to your TORBEAM dev folder> src/OverviewPlot.jl
```