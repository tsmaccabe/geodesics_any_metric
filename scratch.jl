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

(g::Metric)(x) = g.f(x...)

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

function connection(christoffel, du, u)
    F = christoffel(u)
    dims = size(F, 1)
    ddu = zeros(size(du))

    for i in 1:dims for j in 1:dims for k in 1:dims
            ddu[i] -= F[i, j, k] * du[j] * du[k]
    end end end

    return ddu
end


c = 1.
a = 0.3
dims = 2
chart = (u, v) -> [(c + a*cos(u))*cos(v), (c + a*cos(u))*sin(v), a*sin(u)]
torus_metric = Metric((u, v) -> [(c + a*cos(u))^2 0.; 0. a^2])

torus_crstf = (p) -> christoffel(torus_metric, p...)


torus_connection_diffeq = (dq, q, p, t) -> connection(torus_crstf, dq, q)


p = 2π*[rand(), rand()]
x = normalize([rand() for _ in 1:dims])


torus_connection_ode = SecondOrderODEProblem(torus_connection_diffeq, x, p, [0., 20.], saveat=0.01)
torus_connection_sol = solve(torus_connection_ode)

q = [mod.(q_elem[3:4], 2π) for q_elem in torus_connection_sol.u]
q_mat = vcat(q'...)
dq = [q[1:2] for q in torus_connection_sol.u]
dq_mat = vcat(dq'...)

geod = [chart(q_elem[1], q_elem[2]) for q_elem in q]
geod_mat = vcat(geod'...)

scatter(geod_mat[:, 1], geod_mat[:, 2], geod_mat[:, 3], markersize=2., markerstrokewidth=0., xlims=[-a-c, a+c], ylims=[-a-c,a+c], zlims=[-a-c,a+c])