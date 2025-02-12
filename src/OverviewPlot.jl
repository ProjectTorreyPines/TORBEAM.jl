using TORBEAM
using RecipesBase, Plots
using LaTeXStrings
using IMAS
using Interpolations


@recipe function plot_ec_profiles(dd::IMAS.dd)
    seriestype := :path
    xlabel --> L"$rho_\mathrm{tor,norm}$"
    ylabel --> L"$\mathrm{d}P/\mathrm{d}V$ [W M$^{-3}$]"
    legend --> :topright
    for ibeam in 1:length(dd.ec_launchers.beam)
        @series begin
            label := L"\text{Power Beam } $i"
            yaxis := :left
            dd.waves.coherent_wave[ibeam].profiles_1d[].grid.rho_tor_norm, dd.waves.coherent_wave[ibeam].profiles_1d[].power_density
        end
    end
    @series begin
        label := L"\text{Total Power}"
        dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm,
        dd.core_sources.source[source_index].profiles_1d[].electrons.power_inside
         
    end
    @series begin
        ylabel := L"$j$ [A M$^{-2}$]"  # Set the right y-axis label
        nothing, nothing  # Dummy series to attach the label
    end
    for ibeam in 1:length(dd.ec_launchers.beam)
        @series begin
            label := L"\text{Current Beam } $i"
            yaxis := :right
            dd.waves.coherent_wave[ibeam].profiles_1d[].grid.rho_tor_norm, dd.waves.coherent_wave[ibeam].profiles_1d[].current_parallel_density
        end
    end
    
    @series begin
        label := L"\text{Total current}"
        dd.core_sources.source[source_index].profiles_1d[].grid.rho_tor_norm,
        dd.core_sources.source[source_index].profiles_1d[].current_parallel_inside
    end
end


dd = IMAS.json2imas("samples/D3D_170325_trimmed.json")
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.torbeam!(dd, torbeam_params)
if length(dd.waves.coherent_wave) >= 0
    p = plot(plot_ec_profiles(dd))
    savefig("Plots/myplot.pdf")
    display(p)
end