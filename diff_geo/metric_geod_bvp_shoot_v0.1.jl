using ColorSchemes, DifferentialEquations, Interpolations, LinearAlgebra, Plots, SparseArrays, Statistics

include("closest_min.jl")
include("plot_perspetive_v0.2.jl")
include("metric_geodesic_fcn_v0.1.jl")

rot2d(u) = [cos(u) -sin(u); sin(u) cos(u)]

function make_ray(q_fcn, q_ini, a_ini, width, height, t_max=sqrt(2), dt=sqrt(2)/1000, a_max=2pi, da=2pi/15)
    t_rng = (0.:dt:(t_max+dt))
    a_rng = (0.:da:a_max)
    
    ray = [[[0., 0.] for _ in a_rng] for _ in t_rng]
    dray = [[[0., 0.] for _ in a_rng] for _ in t_rng]
    for (ang_i, ang) in enumerate(a_rng)
        w0 = -[cos(ang+a_ini), sin(ang+a_ini)]

        q_ode = SecondOrderODEProblem(q_fcn, w0, q_ini, (0, t_max), saveat=dt)
        q_sol = solve(q_ode, Nystrom4(), dt=dt)
        
        for t_i in axes(t_rng)[1]
            if t_i < length(ray)
                ray[t_i][ang_i] = [mod(q_sol.u[t_i][3], width), mod(q_sol.u[t_i][4], height)]
                dray[t_i][ang_i] = [q_sol.u[t_i][1], q_sol.u[t_i][2]]
            end
        end
    end
    
    return ray[1:end-1], dray[1:end-1], collect(t_rng[1:end-1]), collect(a_rng)
end

c = 1.
a = 0.3
max_dist = sqrt(c^2 + a^2)
t_max = 8*max_dist

dt = 0.01
dang = 0.001*2π

target = 2π.*[rand(), rand()]
r_ini = 2π.*[rand(), rand()]

a_ini = 2π.*rand()

cam_angle = [40, 20]

tol = 0.05

dims = 2
x(u, v) = (c + a*cos(v))*cos(u)
y(u, v) = (c + a*cos(v))*sin(u)
z(v) = a*sin(v)
pullback = (u, v) -> [x(u, v), y(u, v), z(v)]
metric = Metric((u, v) -> [x(u, v)^2 + y(u, v)^2 0.; 0. a^2])

decay_law = (r) -> r^-dims
crstf = (p) -> christoffel(metric, p...)
connection_diffeq = (dq, q, p, t) -> connection(crstf, dq, q, p)

ray, dray, t_rng, a_rng = make_ray(connection_diffeq, r_ini, a_ini, 2pi, 2pi, t_max, dt, 2π, dang)

t_enu = axes(t_rng)[1]
a_enu = axes(a_rng)[1]
dst_grd = zeros(t_enu[end], a_enu[end])
fiber_rotors = zeros(size(a_enu))
for ai in a_enu
    fiber_rotors[ai] = atan(reverse(dray[end][ai])...) - atan(reverse(dray[1][ai])...)
    for ti in t_enu
        dst_grd[ti, ai] = norm(pullback(ray[ti][ai]...) - pullback(target...))
    end
end

min_i = 0
min_idcs = nothing
while isnothing(min_idcs) & (min_i < 100)
    global min_i = min_i + 1
    println(min_i)
    global min_idcs = all_mins(dst_grd, tol)
    if min_i >= 100
        throw("No close approach found."*string(target)*", "*string(r_ini))
    end
end
geod_lengths = [dst_grd[ti, ai] for (ti, ai) in min_idcs]
geod_init_dirs = [dray[1][ai] for (_, ai) in min_idcs]
geod_lie_angles = [fiber_rotors[ai] for (_, ai) in min_idcs]
geod_rotors = rot2d.(geod_lie_angles)

transported = sum(source for source in [
    length^-1 * rotor * fiber for (length, rotor, fiber) in zip(geod_lengths, geod_rotors, geod_init_dirs)
])

for (t_best, ang_best) in min_idcs
    min_θφ = mod.(ray[t_best][ang_best], 2π)

    min_xyz = [x(min_θφ[1], min_θφ[2]), y(min_θφ[1], min_θφ[2]), z(min_θφ[2])]

    p1=heatmap(dst_grd, clims=(minimum(dst_grd), maximum(dst_grd)), color=:vangogh)
    scatter!(p1, [ang_best], [t_best], legend=false, color=:white)

    p2=scatter([r_ini[1], min_θφ[1]], [r_ini[2], min_θφ[2]], xlims=[0.,2π], ylims=[0.,2π], markershape=:x, color=:blue, legend=false)
    scatter!(p2, [target[1]], [target[2]], xlims=[0.,2π], ylims=[0.,2π], markershape=:+, color=:blue, markersize=10)
    ray_plt = zeros((length(t_rng), 2))
    ray_plt = vcat([mod.(ray[ti][ang_best], 2pi) for ti in axes(t_rng)[1]]'...)
    plot!(p2, ray_plt[1:t_best, 1], ray_plt[1:t_best, 2], markersize=1.5, markerstrokewidth=0., markercolor=:orange, marker=:circle, color=:orange)

    ray_3d = hcat([x.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), y.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), z.(ray_plt[1:t_best, 2])]...)

    zi_p_idc = [elem >= 0 for elem in ray_3d[:, 3]]
    zi_p = ray_3d[:, 3][zi_p_idc]

    zi_n_idc = [elem < 0 for elem in ray_3d[:, 3]]
    zi_n = ray_3d[:, 3][zi_n_idc]

    xi_p = ray_3d[:, 1][zi_p_idc]
    xi_n = ray_3d[:, 1][zi_n_idc]
    yi_p = ray_3d[:, 2][zi_p_idc]
    yi_n = ray_3d[:, 2][zi_n_idc]

    p3 = scatter()
    plot_perspective!(p3, c, a, xi_p, yi_p, zi_p, xi_n, yi_n, zi_n, cam_angle=cam_angle)

    ini_xyz = [x(r_ini[1], r_ini[2]), y(r_ini[1], r_ini[2]), z(r_ini[2])]
    targ_xyz = [x(target[1], target[2]), y(target[1], target[2]), z(target[2])]
    scatter!(p3, [targ_xyz[1], ini_xyz[1], min_xyz[1]], [targ_xyz[2], ini_xyz[2], min_xyz[2]], [targ_xyz[3], ini_xyz[3], min_xyz[3]])

    p4 = plot(p1, p2, p3, layout=(3,1), size=(540, 1080))

    display(p4)
    println("Done")
end