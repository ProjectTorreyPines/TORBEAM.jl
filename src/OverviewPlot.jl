using TORBEAM
using RecipesBase, Plots
using LaTeXStrings
using IMAS
using Interpolations


@recipe function plot_ec_profiles(wc1d::IMAS.waves__coherent_wave___profiles_1d; what::Symbol=:power_density)
     @series begin
        wc1d, what
    end
end

@recipe function plot_ec_profiles(wc::IMAS.waves__coherent_wave)
    @series begin
        label := wc.identifier.antenna_name
        wc.profiles_1d[]
    end
end

@recipe function plot_ec_profiles(cwaves::IMAS.IDSvector{<:IMAS.waves__coherent_wave})
    for cwave in cwaves
        @series begin
            cwave
        end
    end
end


@time dd = IMAS.json2imas("samples/D3D_170325_trimmed.json")
dd.global_time = 2.0
torbeam_params = TORBEAM.TorbeamParams()
TORBEAM.torbeam!(dd, torbeam_params)
if length(dd.waves.coherent_wave) >= 0
    p = plot(plot_ec_profiles(dd))
    savefig("Plots/myplot.pdf")
    display(p)
end