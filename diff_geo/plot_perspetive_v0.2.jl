using Plots

include("torus_halves.jl")

cam_angle = [40, 40]

plot_positive!(p, a, b, xi, yi, zi; cam_angle=cam_angle) = scatter!(p, xi, yi, zi, legend=false, camera=cam_angle, xlims=[-(a+2b), (a+2b)], ylims=[-(a+2b), (a+2b)], zlims=[-(a+2b), (a+2b)], aspect_ratio=:equal, markersize=1, markerstrokewidth=0., markercolor=:green)

plot_negative!(p, a, b, xi, yi, zi; cam_angle=cam_angle) = scatter!(p, xi, yi, zi, legend=false, camera=cam_angle, xlims=[-(a+2b), (a+2b)], ylims=[-(a+2b), (a+2b)], zlims=[-(a+2b), (a+2b)], aspect_ratio=:equal, markersize=1, markerstrokewidth=0., markercolor=:blue)

function plot_perspective!(p, a, b, xi_p, yi_p, zi_p, xi_n, yi_n, zi_n; cam_angle=cam_angle)
    elev = cam_angle[2]
    if elev > 0
        torus_half!(p, a, b; top=false)
        plot_negative!(p, a, b, xi_n, yi_n, zi_n; cam_angle=cam_angle)
        torus_half!(p, a, b; top=true)
        plot_positive!(p, a, b, xi_p, yi_p, zi_p; cam_angle=cam_angle)
    else
        torus_half!(p, a, b; top=true)
        plot_positive!(p, a, b, xi_p, yi_p, zi_p; cam_angle=cam_angle)
        torus_half!(p, a, b; top=false)
        plot_negative!(p, a, b, xi_n, yi_n, zi_n; cam_angle=cam_angle)
    end
end