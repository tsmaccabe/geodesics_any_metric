using Colors, DifferentialEquations, Interpolations, LinearAlgebra, Plots, SparseArrays, Statistics

include("closest_min.jl")
include("plot_perspetive.jl")

a = 1.
b = 0.1
torus_max_dist = sqrt(a^2 + b^2)
t_max = 2pi*torus_max_dist

dt = 0.01

dang = 0.01

target = 2π.*[rand(), rand()]

r_ini = [0., 0.]#2π.*[rand(), rand()]
a_ini = 0.#2π.*rand()

cam_angle = [40, 40]

chart_metric(x, y) = sqrt(x^2 + y^2)

x(θ, φ) = (a + b*cos(θ))*cos(φ)
y(θ, φ) = (a + b*cos(θ))*sin(φ)
z(θ) = b*sin(θ)

function torus_metric(u, v, width::Float64=1.0, height::Float64=1.0, n_images::Int=16)
    # Calculate base differences
    Dx = u[1] - v[1]
    Dy = u[2] - v[2]
    
    # Consider n_images in each direction (positive and negative)
    dx_options = [Dx + k*width for k in -n_images:n_images]
    dy_options = [Dy + k*height for k in -n_images:n_images]
    
    # Find minimum distance considering all image combinations
    min_dist = Inf
    for dx in dx_options
        for dy in dy_options
            dist = chart_metric(dx, dy)
            min_dist = min(min_dist, dist)
        end
    end
    
    return min_dist
end

function ddq(w, q, p, t)
    a = p[1]; b = p[2]
    normalize!(w)
    return [-(1/b)*(a + b*cos(q[1]))*sin(q[1])*w[2]^2, 2*b*w[1]*w[2]*sin(q[1])/(a + b*cos(q[1]))]
end

function mod_pos(x, base)
    if x < 0
        quotient = div(abs(x), base)
        x += base*(quotient+1)
    end
    return x % base
end

function make_ray(r_ini, a_ini, width, height, t_max=sqrt(2), dt=sqrt(2)/1000, a_max=2pi, da=2pi/15)
    t_rng = (0.:dt:t_max)
    a_rng = (0.:da:a_max)
    
    ray = [[[0., 0.] for _ in a_rng] for _ in t_rng]
    for (ang_i, ang) in enumerate(a_rng)
        w0 = [cos(ang+a_ini), sin(ang+a_ini)]

        q_ode = SecondOrderODEProblem(ddq, w0, r_ini, (0, t_max), (a, b), saveat=dt)
        q_sol = solve(q_ode, ERKN7())
        
        for (t_i, t) in enumerate(t_rng)
            if t_i < length(ray)
                ray[t_i][ang_i] = q_sol.u[t_i][3:4]
            end
        end
    end
    
    return ray, collect(t_rng), collect(a_rng)
end

ray, t_rng, a_rng = make_ray(r_ini, a_ini, 2pi, 2pi, t_max, dt, 2π, dang)

t_enu = axes(t_rng)[1]
a_enu = axes(a_rng)[1]
dst_grd = zeros(t_enu[end], a_enu[end])
for ti in t_enu
    for ai in a_enu
        dst_grd[ti, ai] = torus_metric(ray[ti][ai], r_ini)
    end
end

min_idc = closest_min(dst_grd, minimum(dst_grd)+5e-2)
t_best, ang_best = min_idc

p1=heatmap(dst_grd, clims=(0, maximum(dst_grd)))
display(p1)