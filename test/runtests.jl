using Plots
using TORBEAM
using TORBEAM.IMAS

using Test


dd = IMAS.hdf2imas("/cscratch/denks/FUSE_dds/170325.h5"; error_on_missing_coordinates=false)
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.run_torbeam(dd, torbeam_params)
p1, p2, p3 = TORBEAM.overview_plot(dd);