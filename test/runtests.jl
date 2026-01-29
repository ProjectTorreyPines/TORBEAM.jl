using Plots
using TORBEAM
using TORBEAM.IMAS

using Test


dd = IMAS.json2imas(dirname(@__DIR__) *"/samples/D3D_170325_trimmed.json"; error_on_missing_coordinates=false)
dd.global_time = 2.0
dd.pulse_schedule.ec.time = zeros(Float64, 1)
dd.pulse_schedule.ec.time[1] = 2.0
resize!(dd.pulse_schedule.ec.beam, length(dd.ec_launchers.beam))
power = 0.0
for i_beam in 1:length(dd.ec_launchers.beam)
    dd.pulse_schedule.ec.beam[i_beam].power_launched.reference = zeros(Float64, 1)
    beam_power = @ddtime(dd.ec_launchers.beam[i_beam].power_launched.data)
    global power += beam_power
    @ddtime(dd.pulse_schedule.ec.beam[i_beam].power_launched.reference = beam_power)
end
dd.pulse_schedule.ec.power_launched.reference = zeros(Float64, 1)
@ddtime(dd.pulse_schedule.ec.power_launched.reference = power)
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.run_torbeam(dd, torbeam_params)
target_power = [569137.0, 505582.0, 520115.0, 
                457320.0, 558655.0, 718337.0]
for i_source in 1:length(dd.core_sources.source)
    @test dd.core_sources.source[i_source].global_quantities[].electrons.power ≈ target_power[i_source] atol=1.0
end