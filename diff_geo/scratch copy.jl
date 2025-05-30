using Interpolations, Plots, SparseArrays, Statistics

chart_metric(x, y) = sqrt(x^2 + y^2)

function torus_metric(u, v, width::Float64=1.0, height::Float64=1.0, n_images::Int=1)
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
        dst_grd[ti, ai] = torus_metric(ray[ti][ai], target)
    end
end

#=dst_itp = linear_interpolation((t_enu, a_enu), dst_grd)

dst_f = (p) -> dst_itp(Int(round(p[1]*(t_enu[end]-1))) + 1, Int(round(p[2]*(a_enu[end]-1))) + 1)=#
function first_blob_center(sparse_matrix)
    m, n = size(sparse_matrix)
    
    # Get row and column indices of non-zero elements
    rows, cols = findnz(sparse_matrix)[1:2]
    
    # Track visited positions to avoid counting same blob twice
    visited = Set{Tuple{Int,Int}}()
    
    # Helper function for BFS
    function explore_blob(start_row, start_col)
        blob_positions = Tuple{Int,Int}[]
        queue = [(start_row, start_col)]
        
        while !isempty(queue)
            row, col = popfirst!(queue)
            
            if (row, col) in visited
                continue
            end
            
            push!(visited, (row, col))
            push!(blob_positions, (row, col))
            
            # Check neighbors (8-connected)
            for dr in -1:1
                for dc in -1:1
                    new_row, new_col = row + dr, col + dc
                    if 1 <= new_row <= m && 1 <= new_col <= n &&
                       sparse_matrix[new_row, new_col] == 1 &&
                       (new_row, new_col) ∉ visited
                        push!(queue, (new_row, new_col))
                    end
                end
            end
        end
        
        return blob_positions
    end
    
    # Search row by row
    for i in 1:m
        for j in 1:n
            if sparse_matrix[i, j] == 1 && (i, j) ∉ visited
                # Found start of a new blob
                blob_positions = explore_blob(i, j)
                
                # Calculate center of blob
                center_row = round(Int, mean([p[1] for p in blob_positions]))
                center_col = round(Int, mean([p[2] for p in blob_positions]))
                
                return (center_row, center_col)
            end
        end
    end
    
    return nothing  # No blob found
end

function closest_min(dst_grd, tol=1e-2)
    dst_mins = sparse(dst_grd .<= tol)
    min_idc = first_blob_center(dst_mins)
    return min_idc
end

#best_val, best_ta = grid_search_2d(dst_f, [t_rng[2], a_ini], [1, 1], [100, 100], 10)

min_idc = closest_min(dst_grd)

using Plots
p1=heatmap(dst_grd, clims=(0., 1.))
scatter!(p1, [min_idc[2]], [min_idc[1]], legend=false, color=:white)
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