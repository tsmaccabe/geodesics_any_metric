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

#dims = 2
#chart = (x, y) -> [x, y]
#euclid_metric = Metric((x, y) -> [x^2 0.; 0. y^2], 2)

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

function connection(christoffel, du, u, p)
    F = christoffel(u)
    dims = size(F, 1)
    ddu = zeros(size(du))

    for i in 1:dims for j in 1:dims for k in 1:dims
            ddu[i] -= F[i, j, k] * du[j] * du[k]
    end end end

    ddu = project_onto_perpendicular(ddu, du)

    return ddu
end


