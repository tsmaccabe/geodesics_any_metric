range = -66:2:66

pitches = vcat(range, reverse(range))

plots = []
anim = @animate for pitch in range
    cam_angle = [40, pitch]
    
    p3 = scatter()
    if pitch < 0
        torus_half!(p3, a, b, top=true)
        torus_half!(p3, a, b)
    else
        torus_half!(p3, a, b)
        torus_half!(p3, a, b, top=true)
    end


    for (t_best, ang_best) in min_idc
        println(cam_angle)
        ray_plt = zeros((length(t_rng), 2))
        ray_plt = vcat([mod_pos.(ray[ti][ang_best], 2pi) for ti in axes(t_rng)[1]]'...)

        ray_3d = hcat([x.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), y.(ray_plt[1:t_best, 1], ray_plt[1:t_best, 2]), z.(ray_plt[1:t_best, 1])]...)

        # Separate indices for z-positive and z-negative parts
        zi_p_idc = [elem >= 0 for elem in ray_3d[:, 3]]
        zi_n_idc = [elem < 0 for elem in ray_3d[:, 3]]
    
        # Collect coordinates for z-negative parts
        append!(all_xi_n, ray_3d[:, 1][zi_n_idc])
        append!(all_yi_n, ray_3d[:, 2][zi_n_idc])
        append!(all_zi_n, ray_3d[:, 3][zi_n_idc])
    
        # Collect coordinates for z-positive parts
        append!(all_xi_p, ray_3d[:, 1][zi_p_idc])
        append!(all_yi_p, ray_3d[:, 2][zi_p_idc])
        append!(all_zi_p, ray_3d[:, 3][zi_p_idc])
    
        min_θφ = mod_pos.(ray[t_best][ang_best], 2π)
        min_xyz = [x(min_θφ[1], min_θφ[2]), y(min_θφ[1], min_θφ[2]), z(min_θφ[1])]

        scatter!(p3, [min_xyz[1]], [min_xyz[2]], [min_xyz[3]], color=:tan, camera=cam_angle)
    end

    if pitch < 0
        scatter!(p3, all_xi_p, all_yi_p, all_zi_p, markersize=0.7, markerstrokewidth=0, color=:green, label="z-positive")
        scatter!(p3, all_xi_n, all_yi_n, all_zi_n, markersize=0.7, markerstrokewidth=0, color=:blue, label="z-negative")
    else
        scatter!(p3, all_xi_n, all_yi_n, all_zi_n, markersize=0.7, markerstrokewidth=0, color=:blue, label="z-negative")
        scatter!(p3, all_xi_p, all_yi_p, all_zi_p, markersize=0.7, markerstrokewidth=0, color=:green, label="z-positive")
    end

    ini_xyz = [x(r_ini[1], r_ini[2]), y(r_ini[1], r_ini[2]), z(r_ini[1])]
    targ_xyz = [x(target[1], target[2]), y(target[1], target[2]), z(target[1])]

    scatter!(p3, [targ_xyz[1], ini_xyz[1]], [targ_xyz[2], ini_xyz[2]], [targ_xyz[3], ini_xyz[3]], color=:tan)
    
    push!(plots, deepcopy(p3))
end

function cleanstring(input_string::String)
    return replace(input_string, r"[[:punct:]]" => "_")
end
using Dates
gif(anim, "torus_geodesics/multi_"*cleanstring(string(now()))*".mp4", fps = 10)
