include("torus_pointtopoint_geodesics_v0.2.jl")
include("torus_halves.jl")

using QuadGK

local_decay(r) = 1/r

function mod_pos(x, base)
    if x < 0
        quotient = div(abs(x), base)
        x += base*(quotient+1)
    end
    return x % base
end

function transport(v::Vector{Float64}, geod_vel_change::Vector{Vector{Float64}})::Vector{Float64}
    # The transported tangent space coordinates are v rotated in the same way as the geodesic's velocity.
    a = atan(geod_vel_change[2], geod_vel_change[1])
    R = [cos(a) -sin(a); sin(a) cos(a)]
    return R*v
end

function gauge_force(geod_fcn::Function, transport::Function, gauge_value::Vector{Float64}, origin::Vector{Float64}, target::Vector{Float64})
    # Unused
    # Allows a gauge field whose fibers are represented by a matrix representation of a symmetry group.

    _, geod_vels, geod_lengths = geod_fcn(origin, target)

    force = [0., 0.]
    for (vel, geod_length) in zip(geod_vels, geod_lengths)
        force += transport.(gauge_value, vel[end]-vel[1]) * local_decay(geod_length)
    end
    return force
end

function gauge_force_radial(geod_vels, geod_lengths)
    force = [0., 0.]
    for (vel, geod_length) in zip(geod_vels, geod_lengths)
        force += vel[end] * local_decay(geod_length)
    end
    return force
end

function ddq(dq, q, p, t)
    null_geod_fcn, force_0 = p
    _, geod_vels, geod_lengths = null_geod_fcn(q[1:3], q[4:6])
    println("t = ", t, ", pos = ", q, ", vels = ", dq)

    ddq1 = [0., force_0*gauge_force_radial(geod_vels, geod_lengths)...]

    ddq2 = [0., force_0*gauge_force_radial([reverse(vel) for vel in geod_vels], [geod_length for geod_length in reverse(geod_lengths)])...]
    println("force = ", [ddq1..., ddq2...])

    return [ddq1..., ddq2...]
end

q10 = [0., 2π*rand(), 2π*rand()]
q20 = [0., 2π*rand(), -2π*rand()]
q0 = [q10..., q20...]

dq10 = [0., 0., 0.]
dq20 = [0., 0., 0.]
dq0 = [dq10..., dq20...]

a = 1.
b = 0.4
a, b = normalize([a, b])
torus_max_dist = sqrt(a^2 + b^2)
s_max = 4pi*torus_max_dist

ds = 0.03#*t_max
dang = 0.003*2π

tol = 0.5
space_manifold_p = (a, b, ds, dang, s_max, tol)

null_geod_fcn = (q1, q2) -> get_null_geods(q1[2:3], q2[2:3], space_manifold_p...)
geods, geod_vels = null_geod_fcn(q10, q20)

geod = geods[1]
geod_vel = geod_vels[1]
t_enu = 1:length(geod)

geod_plt = vcat([mod_pos.(geod[ti], 2pi) for ti in axes(t_enu)[1]]'...)

geod_3d = hcat([x.(geod_plt[:, 1], geod_plt[:, 2]), y.(geod_plt[:, 1], geod_plt[:, 2]), z.(geod_plt[:, 1])]...)

p2 = plot()
torus_half!(p2, a, b)
torus_half!(p2, a, b, top=true)

plot!(p2, geod_3d[:, 1], geod_3d[:, 2], geod_3d[:, 3])

final_vel_3d = hcat([x.(geod_vel[end][1], geod_vel[end][2]), y.(geod_vel[end][1], geod_vel[end][2]), z.(geod_vel[end][1])]...)

plot!(p2, [geod_3d[end, 1]-final_vel_3d[1], geod_3d[end, 1]], [geod_3d[end, 2]-final_vel_3d[2], geod_3d[end, 2]], [geod_3d[end, 3]-final_vel_3d[3], geod_3d[end, 3]], color=:black)