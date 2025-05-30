using ColorSchemes, DifferentialEquations, Interpolations, LinearAlgebra, Plots, SparseArrays, Statistics

include("closest_min.jl")
include("plot_perspetive.jl")

a = 1.
b = 0.8
torus_max_dist = sqrt(a^2 + b^2)
t_max = 5*torus_max_dist

dt = 0.06*t_max

dang = 0.006*2π

# Unstable pairs: 
#   [2.332410484350582, 4.356133081282597], [1.1478757923442422, 2.07371490804144]
#   [5.337753659881493, 4.495935078499196], [3.383490514035348, 6.214457227192285]
target = 2π.*[rand(), rand()]
r_ini = 2π.*[rand(), rand()]

a_ini = 0#2π.*rand()

cam_angle = [40, 20]

tol = 0.1

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
        w0 = -[cos(ang+a_ini), sin(ang+a_ini)]

        q_ode = SecondOrderODEProblem(ddq, w0, r_ini, (0, t_max), (a, b), saveat=dt)
        q_sol = solve(q_ode, SofSpa10(), dt=dt)
        
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
        dst_grd[ti, ai] = torus_metric(ray[ti][ai], target, 2pi, 2pi)
    end
end

min_i = 0
min_idc = nothing
while isnothing(min_idc) & (min_i < 100)
    global min_i = min_i + 1
    println(min_i)
    global min_idc = closest_min(dst_grd, tol)
    if min_i >= 100
        throw("No close approach found."*string(target)*", "*string(r_ini))
    end
end
t_best, ang_best = min_idc

min_θφ = mod_pos.(ray[t_best][ang_best], 2π)
min_xyz = [x(min_θφ[1], min_θφ[2]), y(min_θφ[1], min_θφ[2]), z(min_θφ[1])]

p1=heatmap(dst_grd, clims=(minimum(dst_grd), maximum(dst_grd)), color=:vangogh)
scatter!(p1, [ang_best], [t_best], legend=false, color=:white)

p2=scatter([r_ini[1], min_θφ[1]], [r_ini[2], min_θφ[2]], xlims=[0.,2π], ylims=[0.,2π], markershape=:x, color=:blue, legend=false)
scatter!(p2, [target[1]], [target[2]], xlims=[0.,2π], ylims=[0.,2π], markershape=:+, color=:blue, markersize=10)
ray_plt = zeros((length(t_rng), 2))
ray_plt = vcat([mod_pos.(ray[ti][ang_best], 2pi) for ti in axes(t_rng)[1]]'...)
plot!(p2, ray_plt[1:t_best, 1], ray_plt[1:t_best, 2], markersize=1.5, markerstrokewidth=0., markercolor=:orange, marker=:circle, color=:orange)

ray_3d = hcat([x.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), y.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), z.(ray_plt[1:t_best, 1])]...)

zi_p_idc = [elem >= 0 for elem in ray_3d[:, 3]]
zi_p = ray_3d[:, 3][zi_p_idc]

zi_n_idc = [elem < 0 for elem in ray_3d[:, 3]]
zi_n = ray_3d[:, 3][zi_n_idc]

xi_p = ray_3d[:, 1][zi_p_idc]
xi_n = ray_3d[:, 1][zi_n_idc]
yi_p = ray_3d[:, 2][zi_p_idc]
yi_n = ray_3d[:, 2][zi_n_idc]

p3 = scatter()
plot_perspective!(p3, a, b, xi_p, yi_p, zi_p, xi_n, yi_n, zi_n, cam_angle[2])

ini_xyz = [x(r_ini[1], r_ini[2]), y(r_ini[1], r_ini[2]), z(r_ini[1])]
targ_xyz = [x(target[1], target[2]), y(target[1], target[2]), z(target[1])]
scatter!(p3, [targ_xyz[1], ini_xyz[1], min_xyz[1]], [targ_xyz[2], ini_xyz[2], min_xyz[2]], [targ_xyz[3], ini_xyz[3], min_xyz[3]])

p4 = plot(p1, p2, p3, layout=(3,1), size=(540, 1080))

display(p4)
println("Done")