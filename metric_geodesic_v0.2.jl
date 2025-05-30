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

dims = 2
chart = (x, y) -> [x, y]
euclid_metric = Metric((x, y) -> [x^2 0.; 0. y^2], 2)

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
    #normalize!(du)

    ddu = project_onto_perpendicular(ddu, du)

    return ddu
end

dims = 2
pullback = (x, y) -> [x, y]
metric = Metric((x, y) -> [y^-2 0.; 0. y^-2])

crstf = (p) -> christoffel(metric, p...)

connection_ode_gen = crstf -> ((dq, q, p, t) -> connection(crstf, dq, q))
connection_diffeq = connection_ode_gen(crstf)

p = [rand(), rand()]
x = normalize([rand() for _ in 1:dims])

connection_ode = SecondOrderODEProblem(connection_diffeq, x, p, [0., 5.], saveat=0.03)
connection_sol = solve(connection_ode)

q = [q_elem[3:4] for q_elem in connection_sol.u]
q_mat = vcat(q'...)
dq = [q[1:2] for q in connection_sol.u]
dq_mat = vcat(dq'...)

geod = [pullback(q_elem[1], q_elem[2]) for q_elem in q]
geod_mat = vcat(geod'...)

p1=plot(q_mat[:, 1], q_mat[:, 2])
p2=plot(geod_mat[:, 1], geod_mat[:, 2], mark=:circle, markersize=0.75, markerstrokewidth=0., legend=false)
p3=plot(dq_mat[:, 1], legend=false); plot!(p3, dq_mat[:, 2]); plot!(p3, norm.(dq))
p4=plot(dq_mat[:, 1], dq_mat[:, 2])

p = plot(p1, p2, p3, p4, layout=(4, 1), size=(480, 960))