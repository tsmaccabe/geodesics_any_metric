include("all_mins.jl")
include("plot_perspetive_v0.2.jl")
include("torus_halves.jl")

using ColorSchemes, DifferentialEquations, Interpolations, LinearAlgebra, Plots, SparseArrays, Statistics

chart_metric(x, y) = sqrt(x^2 + y^2)

x(θ, φ) = (a + b*cos(θ))*cos(φ)
y(θ, φ) = (a + b*cos(θ))*sin(φ)
z(θ) = b*sin(θ)

function get_null_geods(origin, target, a, b, dt, dang, t_max, tol=0.01)
    function torus_metric(u, v, width::Float64=1.0, height::Float64=1.0, n_images::Int=16)
        Dx = u[1] - v[1]
        Dy = u[2] - v[2]
        
        dx_options = [Dx + k*width for k in -n_images:n_images]
        dy_options = [Dy + k*height for k in -n_images:n_images]

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
        return [-(1/p[2])*(p[1] + p[2]*cos(q[1]))*sin(q[1])*w[2]^2, 2*p[2]*w[1]*w[2]*sin(q[1])/(p[1] + p[2]*cos(q[1]))]
    end

    function mod_pos(x, base)
        if x < 0
            quotient = div(abs(x), base)
            x += base*(quotient+1)
        end
        return x % base
    end

    function make_ray(origin, a_ini, width, height, t_max=sqrt(2), dt=sqrt(2)/1000, a_max=2pi, da=2pi/15)
        t_rng = (0.:dt:t_max)
        a_rng = (0.:da:a_max)
        
        ray = [[[0., 0.] for _ in t_rng] for _ in a_rng]
        vel = [[[0., 0.] for _ in t_rng] for _ in a_rng]
        for (ang_i, ang) in enumerate(a_rng)
            w0 = -[cos(ang+a_ini), sin(ang+a_ini)]

            q_ode = SecondOrderODEProblem(ddq, w0, origin, (0, t_max), (a, b), saveat=dt)
            q_sol = solve(q_ode, KahanLi8(), dt=dt)
            
            for (t_i, t) in enumerate(t_rng)
                if t_i < length(ray[1])
                    ray[ang_i][t_i] = q_sol.u[t_i][3:4]
                    vel[ang_i][t_i] = q_sol.u[t_i][1:2]
                end
            end
        end
        
        return ray, vel, collect(t_rng), collect(a_rng)
    end

    ray, vel, t_rng, a_rng = make_ray(origin, 0., 2pi, 2pi, t_max, dt, 2π, dang)

    t_enu = axes(t_rng)[1]
    a_enu = axes(a_rng)[1]
    dst_grd = zeros(t_enu[end], a_enu[end])
    for ti in t_enu
        for ai in a_enu
            dst_grd[ti, ai] = torus_metric(ray[ai][ti], target, 2pi, 2pi)
        end
    end

    min_i = 0
    min_idc = nothing
    while isnothing(min_idc) & (min_i < 1)
        min_i = min_i + 1
        min_idc = all_mins(dst_grd, tol)
        if min_i > 1
            throw("No close approach found."*string(target)*", "*string(origin))
        end
    end

    geods = Vector{Vector{Float64}}[]
    tangents = Vector{Vector{Float64}}[]
    geod_lengths = zeros(length(min_idc))
    for (i, (t_best, ang_best)) in enumerate(min_idc)
        push!(geods, ray[ang_best][1:t_best])
        push!(tangents, vel[ang_best][1:t_best])
        geod_lengths[i] = t_best*dt
    end

    return geods, tangents, geod_lengths
end

#=a = 1.
b = 0.4
a, b = normalize([a, b])
torus_max_dist = sqrt(a^2 + b^2)
t_max = 3pi*torus_max_dist

dt = 0.015#*t_max

dang = 0.002*2π

# Unstable pairs: 
#   [2.332410484350582, 4.356133081282597], [1.1478757923442422, 2.07371490804144]
#   [5.337753659881493, 4.495935078499196], [3.383490514035348, 6.214457227192285]
target = 2π.*[rand(), rand()]
origin = 2π.*[rand(), rand()]

a_ini = 0#2π.*rand()

tol = 0.1

geods = get_geods(origin, target, a, b, dt, dang, t_max, tol)
=#