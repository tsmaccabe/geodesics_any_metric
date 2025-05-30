include("torus_pointtopoint_geodesics_v0.2.jl")
include("torus_halves.jl")

using QuadGK

local_decay(r) = 1/r

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

q10 = [0., 2π*0.1, 2π*0.1]
q20 = [0., 2π*0.5, -2π*0.25]
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

tol = 0.02
space_manifold_p = (a, b, ds, dang, s_max, tol)

null_geod_fcn = (q1, q2) -> get_null_geods(q1[2:3], q2[2:3], space_manifold_p...)
geods, geod_vels = null_geod_fcn(q10, q20)

force_0 = 1.
p = (null_geod_fcn, force_0)

dl = 0.1 # Proper Time
l_max = 5.

q_ode = SecondOrderODEProblem(ddq, dq0, q0, (0, l_max), p, saveat=dl)
q_sol = solve(q_ode, DPRKN8())
q = [q_sol.u[si][7:12] for si in axes(q_sol.u)[1]]

q1 = vcat([collect(q[si][1:3]') for si in axes(q)[1]]...)
q2 = vcat([collect(q[si][4:6]') for si in axes(q)[1]]...)

p1 = plot(q1[:, 2], q1[:, 3], legend=false)
plot!(p1, q2[:, 2], q2[:, 3])
scatter!(p1, [q1[end, 2], q2[end, 2]], [q1[end, 3], q2[end, 3]])
display(p1)

q1_space_embedded = hcat(x.(q1[:, 2], q1[:, 3]), y.(q1[:, 2], q1[:, 3]), z.(q1[:, 3]))
q2_space_embedded = hcat(x.(q2[:, 2], q2[:, 3]), y.(q2[:, 2], q2[:, 3]), z.(q2[:, 3]))

p2 = plot()
torus_half!(p2, a, b)
torus_half!(p2, a, b, top=true)
plot!(p2, q1_space_embedded[:, 1], q1_space_embedded[:, 2], q1_space_embedded[:, 3], linewidth=5.)
plot!(p2, q2_space_embedded[:, 1], q2_space_embedded[:, 2], q2_space_embedded[:, 3], linewidth=5.)
scatter!(p2, [q1_space_embedded[end, 1]], [q1_space_embedded[end, 2]], [q1_space_embedded[end, 3]], markersize=5.)
scatter!(p2, [q2_space_embedded[end, 1]], [q2_space_embedded[end, 2]], [q2_space_embedded[end, 3]], markersize=5.)
display(p2)