using TORBEAM
using Plots
using IMAS
using Interpolations

dd = IMAS.json2imas("samples/D3D_170325_trimmed.json")
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.torbeam!(dd, torbeam_params)
if length(dd.waves.coherent_wave) >= 0
    p = plot( layout=(2, 1))
    total_power = zeros(Float64, length(dd.waves.coherent_wave[1].profiles_1d[1].power_density))
    total_current = zeros(Float64, length(dd.waves.coherent_wave[1].profiles_1d[1].current_parallel_density))
    axis = dd.waves.coherent_wave[1].profiles_1d[1].grid.rho_tor_norm
    for ibeam =1:length(dd.ec_launchers.beam)
        plot!(p, dd.waves.coherent_wave[ibeam].profiles_1d[].grid.rho_tor_norm, 
              [dd.waves.coherent_wave[ibeam].profiles_1d[].power_density, 
               dd.waves.coherent_wave[ibeam].profiles_1d[].current_parallel_density], layout=(2, 1))
    end
    source_index = 1
    plot!(p, dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm, 
          [dd.core_sources.source[source_index].profiles_1d[].electrons.power_inside, 
           dd.core_sources.source[source_index].profiles_1d[].current_parallel_inside], layout=(2, 1))
    savefig("Plots/myplot.pdf")
    display(p)
end