module TORBEAMext

using TORBEAM
using Plots

function overview_plot(dd)
    p_rz = plot()
    p_xy = plot()
    p_prof = plot(; layout=(2, 1))
    plot!(p_rz, dd.equilibrium.time_slice[]; cx=true)

    phi = LinRange(0.0, 2.0 * pi, 200)
    ipsi = findmin(abs.(dd.equilibrium.time_slice[].profiles_1d.psi .- 1.0))[2]
    x_inboard = dd.equilibrium.time_slice[].profiles_1d.r_inboard[ipsi] .* cos.(phi)
    y_inboard = dd.equilibrium.time_slice[].profiles_1d.r_inboard[ipsi] .* sin.(phi)
    x_outboard = dd.equilibrium.time_slice[].profiles_1d.r_outboard[ipsi] .* cos.(phi)
    y_outboard = dd.equilibrium.time_slice[].profiles_1d.r_outboard[ipsi] .* sin.(phi)
    plot!(p_xy, x_inboard, y_inboard; cx=true)
    plot!(p_xy, x_outboard, y_outboard; cx=true)

    for ibeam in 1:length(dd.waves.coherent_wave)
        plot!(p_rz, dd.waves.coherent_wave[ibeam])
        plot!(p_xy, dd.waves.coherent_wave[ibeam]; top=true)
        plot!(p_prof, dd.waves.coherent_wave[ibeam].profiles_1d; xlimits=(0.1, 0.4))
    end
    savefig(p_rz, "Plots/TORBEAM_rz.pdf")
    savefig(p_xy, "Plots/TORBEAM_xy.pdf")
    savefig(p_prof, "Plots/TORBEAM_profiles.pdf")
    return p_rz, p_xy, p_prof
end

end  # module TORBEAMext
