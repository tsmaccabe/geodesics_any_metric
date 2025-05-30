import AbstractDifferentiation as AD, Zygote
backend = AD.ForwardDiffBackend()
jacobian = (metric, xi) -> AD.jacobian(backend, metric, xi)[1]

using DifferentialEquations, LinearAlgebra, Plots

mutable struct Metric <: Function
    dims::Integer
    f::Function
    function Metric(f, dims=2)
        return new(dims, f)
    end
end

(g::Metric)(x...) = g.f(x...)
(g::Metric)(x::Vector) = g.f(x...)

mutable struct Manifold
    dims::Integer
    time::Bool
    metric::Metric
    function Manifold(metric, time=false)
        return new(size(metric(0., 0.)), time, metric) 
    end
end

jacob = (metric::Metric, x, y) -> reshape(jacobian(metric, [x, y]), (2*ones(Integer, 3))...)

function christoffel(metric::Metric, x, y)
    Ginv = metric(x, y)^-1
    J = jacob(metric, x, y)
    F = zeros(typeof(x), (dims*ones(Integer, 3))...)

    for i in 1:dims for j in 1:dims for k in 1:dims for l in 1:dims
            F[i, j, k] += 0.5 * Ginv[i, l] * (J[l, j, k] + J[l, k, j] - J[j, k, l])
    end end end end

    return F
end

function project_onto_perpendicular(ddu, du)
    projection = (dot(ddu, du) / dot(du, du)) * du
    
    perpendicular_component = ddu - projection
    
    return perpendicular_component
end

function connection(christoffel, du, u)
    F = christoffel(u)
    dims = size(F, 1)
    ddu = zeros(size(du))

    for i in 1:dims for j in 1:dims for k in 1:dims
            ddu[i] -= F[i, j, k] * du[j] * du[k]
    end end end

    normalize!(du)
    ddu = project_onto_perpendicular(ddu, du)

    return ddu
end

c = 1.
a = 0.8
dims = 2
pullback = (u, v) -> [(c + a*cos(v))*cos(u), (c + a*cos(v))*sin(u), a*sin(v)]
metric = Metric((u, v) -> [(c + a*cos(v))^2 0.; 0. a^2])

crstf = (p) -> christoffel(metric, p...)


connection_diffeq = (dq, q, p, t) -> connection(crstf, dq, q)


q0 = 2π*[rand(), rand()]
q2 = 2π*[rand(), rand()]
x = normalize([rand() for _ in 1:dims])

connection_ode = SecondOrderODEProblem(connection_diffeq, x, q0, [0., 50.], saveat=0.03)
connection_sol = solve(connection_ode)

q = [mod.(q_elem[3:4], 2π) for q_elem in connection_sol.u]
q_mat = vcat(q'...)
dq = [q[1:2] for q in connection_sol.u]
dq_mat = vcat(dq'...)

geod = [pullback(q_elem[1], q_elem[2]) for q_elem in q]
geod_mat = vcat(geod'...)


function bvp_shoot(connection_diffeq, q1, q2, shots::Integer, max_dist, layers::Integer=3, metric=Metric((x, y) -> [1. 0.; 0. 1.], 2), tol=1e-1, saveat=5e-2)
    dθ = 2π/shots
    θ = 0.:dθ:(2π-dθ)

    chart_distance = (q1) -> sqrt(det(metric(q1...)))

    shortest_geod = nothing

    for l in 1:layers
        min_dists = Float64[]
        geods = Matrix{Float64}[]

        for θ_i in θ
            dq1 = -[cos(θ_i), sin(θ_i)]

            connection_ode = SecondOrderODEProblem(connection_diffeq, dq1, q1, [0., max_dist], saveat=0.03)
            connection_sol = solve(connection_ode)

            dispq = [mod.(q_elem[3:4], 2π) .- q2 for q_elem in connection_sol.u]
            disp_mat = vcat(dispq'...)
            push!(geods, disp_mat)

            dists = chart_distance.(dispq)
            push!(min_dists, minimum(dists))
        end
        
        shortest_i = argmin(min_dists)

        if min_dists[shortest_i] < tol
            shortest_geod = geods[shortest_i]
            break
        end

        shortest_θ = θ[shortest_i]
        if l == 1
            next_θ = θ[mod(shortest_i+1, shots)+1]
            prev_θ = θ[mod(shortest_i-3, shots)+1]
        else
            next_θ = θ[minimum([shortest_i+2, shots])]
            prev_θ = θ[maximum([shortest_i-2, 1])]
        end
        range_θ = next_θ - prev_θ
        println(min_dists[shortest_i])
        
        dθ = range_θ/shots
        θ = prev_θ:dθ:next_θ
    end

    return shortest_geod
end

qp1 = pullback(q0...)
qp2 = pullback(q2...)

geod = bvp_shoot(connection_diffeq, q2, q0, 1000, 1., 1, metric, 0.05)
geod = [pullback(q_elem[1], q_elem[2]) for q_elem in q]
geod_mat = vcat(geod'...)

p1=plot(q_mat[:, 1], mod.(q_mat[:, 2] .+ π, 2π))
p2=plot(geod_mat[:, 1], geod_mat[:, 2], geod_mat[:, 3], mark=:circle, markersize=0.75, markerstrokewidth=0., xlims=[-a-c, a+c], ylims=[-a-c,a+c], zlims=[-a-c,a+c], legend=false)
scatter!(p2, [qp1[1], qp2[1]], [qp1[2], qp2[2]], [qp1[3], qp2[3]], markersize=2.)
p3=plot(dq_mat[:, 1], legend=false); plot!(p3, dq_mat[:, 2]); plot!(p3, norm.(dq))

p = plot(p1, p2, p3, layout=(3, 1), size=(480, 960))
display(p)