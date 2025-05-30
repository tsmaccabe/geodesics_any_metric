using Dates, DifferentialEquations, LinearAlgebra, Plots

a = 10.
b = 9.9

h0 = pi/2
f0 = 0.
q0 = [h0, f0]

vh0 = 1.
vf0 = 1. /(a-b)
w0 = normalize([vh0, vf0])

function ddq(w, q, p, t)
    a = p[1]; b = p[2]
    return [-(1/b)*(a + b*cos(q[1]))*sin(q[1])*w[2]^2, 2*b*w[1]*w[2]*sin(q[1])/(a + b*cos(q[1]))]
end

dt = 0.01
q_ode = SecondOrderODEProblem(ddq, w0, q0, (0, 1000), (a, b), saveat=dt)
q_sol = solve(q_ode, SofSpa10(), dt=dt)

t = q_sol.t
hi = [ui[3] for ui in q_sol.u]
fi = [ui[4] for ui in q_sol.u]

x(θ, φ) = (a + b*cos(θ))*cos(φ)
y(θ, φ) = (a + b*cos(θ))*sin(φ)
z(θ) = b*sin(θ)

θ(z) = asin(z/b)
φ(x, z) = acos(clamp(x/(a + b*cos(asin(z/b))), -1, 1))

xi = x.(hi, fi)
yi = y.(hi, fi)
zi = z.(hi)

zi_p_idc = [elem >= 0 for elem in zi]
zi_p = zi[zi_p_idc]

zi_n_idc = [elem < 0 for elem in zi]
zi_n = zi[zi_n_idc]

xi_p = xi[zi_p_idc]
xi_n = xi[zi_n_idc]
yi_p = yi[zi_p_idc]
yi_n = yi[zi_n_idc]

θi = θ.(zi)
φi = φ.(xi, zi)

plot(θi, φi)

cam_angle = [40, 40]

plot_positive!(p) = scatter!(p, xi_n, yi_n, zi_n, legend=false, camera=cam_angle, xlims=[-(a+b), (a+b)], ylims=[-(a+b), (a+b)], zlims=[-(a+b), (a+b)], aspect_ratio=:equal, markersize=0.6, markerstrokewidth=0., markercolor=:blue)

plot_negative!(p) = scatter!(p, xi_p, yi_p, zi_p, legend=false, camera=cam_angle, xlims=[-(a+b), (a+b)], ylims=[-(a+b), (a+b)], zlims=[-(a+b), (a+b)], aspect_ratio=:equal, markersize=0.6, markerstrokewidth=0., markercolor=:orange)

function plot_perspective!(p, elev)
    if elev > 0
        plot_positive!(p)
        plot_negative!(p)
    else
        plot_negative!(p)
        plot_positive!(p)
    end
end

p = scatter()
plot_perspective!(p, cam_angle[2])
#scatter!(p, [xi[1]], [yi[1]], [zi[3]], markercolor=:white, markersize=4)
display(p)

function cleanstring(input_string::String)
    return replace(input_string, r"[[:punct:]]" => "_")
end
savefig("field_automaton_0.2\\geodesic_images\\torus_geodesic"*cleanstring(string(now()))*".png")