using TORBEAM
using TORBEAM.IMAS
using Test
using Plots

dd = IMAS.json2imas(joinpath(@__DIR__, "..", "samples", "D3D_170325_trimmed.json"))
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.run_torbeam(dd, torbeam_params)
p1, p2, p3 = TORBEAM.overview_plot(dd);