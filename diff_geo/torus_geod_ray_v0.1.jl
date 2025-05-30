using Dates, DifferentialEquations, Plots

a = 10.
b = 5

f0 = 0
h0 = 0

vf0 = 1
vh0 = -1/2

vf(t) = vf0
vh(h, p, t) = vh0 + vf0*log(p[1] + p[2]*cos(h))


vh_ode = ODEProblem(vh, (vh0, h0), (0, 1), (a, b), saveat=0.01)
vh_sol = solve(vh_ode)

mod_sym(num, base) = (num + base*(num < 0)) % base

t = vh_sol.t
fi = mod_sym.(vf0*t, 2π)
hi = mod_sym.(vh_sol.u, 2π)

x(θ, φ) = (a + b*cos(θ))*cos(φ)
y(θ, φ) = (a + b*cos(θ))*sin(φ)
z(θ) = b*sin(θ)

xi = x.(hi, fi)
yi = y.(hi, fi)
zi = z.(hi)

zi_p = [elem*(elem > 0) for elem in zi]
zi_p_idc = (zi_p .!= 0)
zi_p = zi_p[zi_p_idc]

zi_n = [elem*(elem <= 0) for elem in zi]
zi_n_idc = (zi_n .!= 0)
zi_n = zi_n[zi_n_idc]

xi_p = xi[zi_p_idc]
xi_n = xi[zi_n_idc]
yi_p = yi[zi_p_idc]
yi_n = yi[zi_n_idc]

cam_angle = [40, 40]

plot_positive!(p) = scatter!(p, xi_n, yi_n, zi_n, legend=false, camera=cam_angle, xlims=[-15, 15], ylims=[-15, 15], zlims=[-15, 15], aspect_ratio=:equal, markersize=0.6, markerstrokewidth=0., markercolor=:blue)

plot_negative!(p) = scatter!(p, xi_p, yi_p, zi_p, legend=false, camera=cam_angle, xlims=[-15, 15], ylims=[-15, 15], zlims=[-15, 15], aspect_ratio=:equal, markersize=0.6, markerstrokewidth=0., markercolor=:orange)

function plot_perspective(elev)
    p=scatter()

    if elev > 0
        plot_positive!(p)
        plot_negative!(p)
    else
        plot_negative!(p)
        plot_positive!(p)
    end

    return p
end

p = plot_perspective(cam_angle[2])
scatter!(p, [xi[1]], [yi[1]], [zi[3]], markersize=2)
display(p)

function cleanstring(input_string::String)
    return replace(input_string, r"[[:punct:]]" => "_")
end
savefig("field_automaton_0.2\\geodesic_images\\torus_geodesic"*cleanstring(string(now()))*".png")