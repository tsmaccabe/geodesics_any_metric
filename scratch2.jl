using LinearAlgebra, DifferentialEquations, Plots

u0_step = 2π/100
u0 = collect(0:u0_step:(2π-u0_step))

r0 = [0., 0.]

Q = [r0 for _ in u0]
dQ = copy(u0)

function distance(r1, r2)
    return norm(r1 - r2)
end

function split_condition(Q, threshold=0.1)
    split_ids = []
    split_dists = []
    for (i, (q_this, q_next)) in enumerate(zip(Q, vcat(Q[2:end], Q[1])))
        dist = distance(q_next, q_this)
        if dist > threshold
            push!(split_ids, i)
            push!(split_dists, dist)
        end
    end
    return split_ids, split_dists
end

function split_trajectories(Q, threshold=0.1)
    split_ids = split_condition(Q, threshold)
    r0_splits = []
    for i in reverse(split_ids)
        insert!(Q, i + 1, (Q[i] + circshift(Q, 1)[i]) / 2)
    end
    return r0_splits
end