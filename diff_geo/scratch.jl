using Interpolations

chart_metric(x, y) = sqrt(x^2 + y^2)

function chart_dst(u, v, width::Float64=1.0, height::Float64=1.0)
    Dx = abs(u[1] - v[1])
    Dy = abs(u[1] - v[2])
    
    Dx = min(Dx, width - Dx)
    Dy = min(Dy, height - Dy)
    
    return chart_metric(Dx, Dy)
end

target = [rand(), rand()]

r_ini = [rand(), rand()]
a_ini = rand()*2pi

function make_ray(r_ini, a_ini, width, height, t_max=sqrt(2), dt=sqrt(2)/100)
    t_rng = (0:dt:t_max)
    ray_points = hcat((r_ini[1] .+ cos(a_ini)*collect(t_rng)) .% width, (r_ini[2] .+ sin(a_ini)*collect(t_rng)) .% height)

    itpx = linear_interpolation(t_rng, [x for x in ray_points[:, 1]])
    itpy = linear_interpolation(t_rng, [y for y in ray_points[:, 2]])

    return (t) -> [itpx(t), itpy(t)], t_rng
end
ray, t_rng = make_ray(r_ini, a_ini, 1, 1)

to_target(t) = chart_dst(target, ray(t))

function grid_search_1d(f::Function, x_ini::Real, x_fin::Real, x_num::Integer, scale_iter_num::Integer)::Real
    best_x = x_ini
    best_val = f(x_ini)

    x_ini_p = copy(x_ini)
    x_fin_p = copy(x_fin)
    
    for scale_iter in 1:scale_iter_num
        dx = (x_fin_p - x_ini_p) / x_num
        x_points = collect(x_ini_p:dx:x_fin_p)
        
        # Store all points and their values
        points_vals = [(x, f(x)) for x in x_points]
        
        # Sort by function value to find 3 best points
        sort!(points_vals, by=x->x[2])
        three_best = points_vals[1:min(3, length(points_vals))]
        println(points_vals)
        println(three_best[1][2])
        
        # Update global best
        if three_best[1][2] < best_val
            best_val = three_best[1][2]
            best_x = three_best[1][1]
        end
        # Find range covering 3 best points
        min_x = minimum(p[1] for p in three_best)
        max_x = maximum(p[1] for p in three_best)
        
        # Set new search range around these points with some padding
        padding = (max_x - min_x) * 0.2  # 20% padding
        x_ini_p = max(min_x - padding, x_ini_p)
        x_fin_p = min(max_x + padding, x_fin_p)
    end
    
    return best_val
end

x_ini = 0.#chart_dst(target, r_ini)/2
println(grid_search_1d(to_target, x_ini, sqrt(2), 20, 5))

using Plots
p=scatter([r_ini[1], target[1]], [r_ini[2], target[2]], xlims=[0.,1.], ylims=[0.,1.])
ray_plt = ray.(t_rng)
ray_plt = hcat([elem[1] for elem in ray_plt], [elem[2] for elem in ray_plt])
plot!(p, ray_plt[:, 1], ray_plt[:, 2])
#=a_min_guess = 0.

t_tol = 1e-3
dst_prv = Inf=#
