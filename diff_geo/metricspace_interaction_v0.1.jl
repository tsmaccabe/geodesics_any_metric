using ColorSchemes, DifferentialEquations, Interpolations, LinearAlgebra, Plots, SparseArrays, Statistics

include("closest_min.jl")
include("plot_perspetive_v0.2.jl")
include("metric_geodesic_fcn_v0.1.jl")

rot2d(u) = [cos(u) -sin(u); sin(u) cos(u)]

function make_ray(q_fcn, q_ini, a_ini, t_max=sqrt(2), dt=sqrt(2)/1000, a_max=2π, da=2pi/15, width=2π, height=2π)
    t_rng = (0.:dt:(t_max+dt))
    a_rng = (0.:da:a_max)
    
    ray = [[[0., 0.] for _ in a_rng] for _ in t_rng]
    dray = [[[0., 0.] for _ in a_rng] for _ in t_rng]
    for (ang_i, ang) in enumerate(a_rng)
        if ang_i % 10 == 0
            println("ray ", ang_i)
        end
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

function lie_multi_transporter(start, target, metric, pullback, t_max, dt, tol=0.3, angle_init=0.; decay_law=(r)->r^-1, include_neg_transports=true)
    crstf = (u) -> christoffel(metric, u...)
    connection_ode = (dq, q, p, t) -> connection(crstf, dq, q, p)

    ray, dray, t_rng, a_rng = make_ray(connection_ode, start, angle_init, t_max, dt, 2π, dang)

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

    #min_i = 0
    min_idcs = nothing
    #while isnothing(min_idcs) & (min_i < 100)
        #min_i = min_i + 1
        #println(min_i)
        min_idcs = all_mins(dst_grd, tol)
        println(min_idcs)
        if length(min_idcs) == 0
            throw("No close approach found."*string(target)*", "*string(start))
        end
    #end
    geod_lengths = [dst_grd[ti, ai] for (ti, ai) in min_idcs]
    geod_init_dirs = [dray[1][ai] for (_, ai) in min_idcs]
    geod_lie_angles = [fiber_rotors[ai] for (_, ai) in min_idcs]
    geod_rotors = rot2d.(geod_lie_angles)

    transported_vec = [decay_law(length) * rotor * init_dir for (length, rotor, init_dir) in zip(geod_lengths, geod_rotors, geod_init_dirs)]
    if length(min_idcs) > 0
        transported = gauge_scale*sum(transported_vec)
    else
        transported = [0., 0.] # No connections found (for smooth manifolds, this is probably a sign of an error)
    end
    
    geod_last_dirs = [dray[end][ai] for (_, ai) in min_idcs]
    geod_rotors_neg = rot2d.(-geod_lie_angles)

    transported_neg_vec = [decay_law(length) * rotor * -last_dir for (length, rotor, last_dir) in zip(geod_lengths, geod_rotors_neg, geod_last_dirs)]
    if length(min_idcs) > 0
        transported_neg = gauge_scale*sum(transported_neg_vec)
    else
        transported_neg = [0., 0.] # No connections found (for smooth manifolds, this is probably a sign of an error)
    end

    return transported, transported_neg
end

c = 1.
a = 0.4
max_dist = sqrt(c^2 + a^2)

a_ini = 2π.*rand()
tol = 0.1

dims = 2
x(u, v) = (c + a*cos(v))*cos(u)
y(u, v) = (c + a*cos(v))*sin(u)
z(v) = a*sin(v)
pullback = (u, v) -> [x(u, v), y(u, v), z(v)]
metric = Metric((u, v) -> [x(u, v)^2 + y(u, v)^2 0.; 0. a^2])

t_max = 10.

dt_gauge = 0.1
dang = 0.05*2π

gauge_scale = 0.1 

function ddq(dq, q, p, t)
    # p should be (gauge_scale, metric, pullback, t_max, dt_gauge, tol, a_ini)
    println(dq, ", ", q, ", ", t)
    transported, transported_neg = lie_multi_transporter(q[1:2], q[3:4], p[2:7]...)
    force = p[1]*[transported..., transported_neg...]
    println(force)
    return force
end

dt_sim = 1e-2
dt_sim_save = dt_sim
sim_time = 10.
x0 = 2π.*[0.0, -0.1, 0.1, 0.1]
dx0 = zeros(4)

particles_ode = SecondOrderODEProblem(ddq, dx0, x0, (0., sim_time), (gauge_scale, metric, pullback, t_max, dt_gauge, tol, a_ini), saveat=dt_sim_save)
particles_sol = solve(particles_ode, dt=dt_sim)

q1 = [u[5:6] for u in particles_sol.u]
q2 = [u[7:8] for u in particles_sol.u]

q1_mat = hcat(q1...)
q2_mat = hcat(q2...)

p1=plot()
plot!(p1, q1_mat[1, :], q1_mat[2, :])
plot!(p1, q2_mat[1, :], q2_mat[2, :])
display(p1)

r1 = [pullback(q...) for q in q1]
r2 = [pullback(q...) for q in q2]

r1_mat = hcat(r1...)
r2_mat = hcat(r2...)

p2=plot()
plot!(p2, r1_mat[1, :], r1_mat[2, :], r1_mat[3, :])
plot!(p2, r2_mat[1, :], r2_mat[2, :], r2_mat[3, :])
display(p2)