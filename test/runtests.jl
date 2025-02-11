using TORBEAM
using Test
using IMAS

dd = IMAS.json2imas("samples/D3D_170325_trimmed.json")
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.torbeam!(dd, torbeam_params)