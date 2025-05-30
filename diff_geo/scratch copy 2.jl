using Interpolations, Plots

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

function make_ray(r_ini, a_ini, width, height, t_max=sqrt(2), dt=sqrt(2)/1000, a_max=2pi, da=2pi/1000)
    t_rng = (0.:dt:t_max)
    a_rng = (a_ini:da:(a_max+a_ini))

    ray = [[[(r_ini[1] + cos(a)*t) .% width, (r_ini[2] + sin(a)*t) .% height] for a in a_rng] for t in t_rng]
    
    return ray, collect(t_rng), collect(a_rng)
end
ray, t_rng, a_rng = make_ray(r_ini, a_ini, 1, 1)

t_enu = axes(t_rng)[1]
a_enu = axes(a_rng)[1]
dst_grd = zeros(t_enu[end], a_enu[end])
for ti in t_enu
    for ai in a_enu
        dst_grd[ti, ai] = chart_dst(ray[ti][ai], target)
    end
end

dst_itp = linear_interpolation((t_rng, a_rng), dst_grd)

dst_f = (p) -> dst_itp(p[1], p[2])

function grid_search_2d(f, x_ini::Vector{<:Real}, x_fin::Vector{<:Real}, 
    x_num::Vector{<:Integer}, scale_iter_num::Integer)::Real
    @assert length(x_ini) == 2 && length(x_fin) == 2 && length(x_num) == 2 "Input dimensions must be 2D"

    best_x = copy(x_ini)
    best_val = f(x_ini)

    x_ini_p = copy(x_ini)
    x_fin_p = copy(x_fin)

    for scale_iter in 1:scale_iter_num
        # Calculate step sizes for both dimensions
        dx = (x_fin_p - x_ini_p) ./ x_num

        # Create grid points
        x1_points = collect(x_ini_p[1]:dx[1]:x_fin_p[1])
        x2_points = collect(x_ini_p[2]:dx[2]:x_fin_p[2])

        # Store all points and their values
        points_vals = [(Float64[x1, x2], f(x1, x2)) 
            for x1 in x1_points, x2 in x2_points]
                points_vals = reshape(points_vals, :)  # Flatten to 1D array

                # Sort by function value to find 3 best points
                sort!(points_vals, by=x->x[2])
                three_best = points_vals[1:min(3, length(points_vals))]

                println("Iteration $scale_iter best value: $(three_best[1][2])")

                # Update global best
                if three_best[1][2] < best_val
                    best_val = three_best[1][2]
                    best_x = three_best[1][1]
                end
                println(best_x)

                # Find ranges covering 3 best points for each dimension
                min_x = minimum(p[1][1] for p in three_best)
                max_x = maximum(p[1][1] for p in three_best)
                min_y = minimum(p[1][2] for p in three_best)
                max_y = maximum(p[1][2] for p in three_best)

                # Set new search ranges around these points with padding
                padding_x = (max_x - min_x) * 0.2
                padding_y = (max_y - min_y) * 0.2

                x_ini_p = [
                max(min_x - padding_x, x_ini[1]),
                max(min_y - padding_y, x_ini[2])
                ]
                x_fin_p = [
                min(max_x + padding_x, x_fin[1]),
                min(max_y + padding_y, x_fin[2])
        ]
    end

    return best_val
end

#best_val, best_ta = grid_search_2d(dst_f, [t_rng[2], a_ini], [1, 1], [100, 100], 10)



using Plots
p1=heatmap(dst_grd, clims=(0., 1.))
display(p1)
p2=scatter([r_ini[1], target[1]], [r_ini[2], target[2]], xlims=[0.,1.], ylims=[0.,1.])
plot!(p2, ray_plt[:, 1], ray_plt[:, 2])
# Create a single figure with vertical layout
plot(p1, p2, layout=(2,1))
#=a_min_guess = 0.

t_tol = 1e-3
dst_prv = Inf=#

x_ini = 0.#chart_dst(target, r_ini)/2
println(grid_search_2d(dst_f, [x_ini, a_ini], [sqrt(2), 2pi], [20, 20], 5))