using ColorSchemes, Plots

function torus_half!(p, c, a, Δα=0.2, θ_dmn=0:Δα/3:2π, φ_dmn=0:Δα/2:2π; top=false)
    x(u, v) = (c + a*cos(v))*cos(u)
    y(u, v) = (c + a*cos(v))*sin(u)
    z(v) = a*sin(v)

    θφ_rng = vec(collect(Iterators.product(θ_dmn, φ_dmn)))
    θ_rng = [θφ[1] for θφ in θφ_rng]
    φ_rng = [θφ[2] for θφ in θφ_rng]

    r_rng = sort([(x(θ, φ), y(θ, φ), z(φ)) for (θ, φ) in zip(θ_rng, φ_rng)], by = r -> r[3])
    x_rng = [r[1] for r in r_rng]
    y_rng = [r[2] for r in r_rng]
    z_rng = [r[3] for r in r_rng]
    if top
        x_rng = x_rng[z_rng .> 0]
        y_rng = y_rng[z_rng .> 0]
        z_rng = z_rng[z_rng .> 0]
    else
        x_rng = x_rng[z_rng .< 0]
        y_rng = y_rng[z_rng .< 0]
        z_rng = z_rng[z_rng .< 0]
    end

    ε = 0 # stability offset
    R_rng = (x_rng .- minimum(x_rng)) / (maximum(x_rng) - minimum(x_rng))
    G_rng = (y_rng .- minimum(y_rng)) / (maximum(y_rng) - minimum(y_rng))
    B_rng = (z_rng .- minimum(z_rng)) / (maximum(z_rng) - minimum(z_rng))
    RGB_rng = [RGB(R, a, G)/255 for (R, G, B) in zip(R_rng, G_rng, B_rng)]

    if top
        color = :brown
    else
        color = :grey
    end
    scatter!(p, x_rng, y_rng, z_rng, color = color, markersize = 1, markerstrokewidth = 0., camera = (30, 30), xlims = [-2, 2], ylims = [-2, 2], zlims = [-2, 2], legend=false)
end