using Plots

include("torus_halves.jl")

cam_angle = [40, 40]

plot_positive!(p, xi, yi, zi) = scatter!(p, xi, yi, zi, legend=false, camera=cam_angle, xlims=[-(a+b), (a+b)], ylims=[-(a+b), (a+b)], zlims=[-(a+b), (a+b)], aspect_ratio=:equal, markersize=1, markerstrokewidth=0., markercolor=:orange)

plot_negative!(p, xi, yi, zi) = scatter!(p, xi, yi, zi, legend=false, camera=cam_angle, xlims=[-(a+b), (a+b)], ylims=[-(a+b), (a+b)], zlims=[-(a+b), (a+b)], aspect_ratio=:equal, markersize=1, markerstrokewidth=0., markercolor=:blue)

function plot_perspective!(p, a, b, xi_p, yi_p, zi_p, xi_n, yi_n, zi_n, elev)
    if elev > 0
        plot_negative!(p, xi_n, yi_n, zi_n)
        plot_positive!(p, xi_p, yi_p, zi_p)
    else
        plot_positive!(p, xi_p, yi_p, zi_p)
        plot_negative!(p, xi_n, yi_n, zi_n)
    end
end