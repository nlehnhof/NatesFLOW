using Plots
using LinearAlgebra

# Compute induced velocity due to vortex
function induced_velocity(gamma, r)
    return cross(gamma, r) / (2π * norm(r)^2)
end

# Compute r vector
function compute_r(point1, point2)
    x1, y1 = point1
    x2, y2 = point2
    return [x2 - x1, y2 - y1, 0.0]
end

# Compute total velocity at a point due to all vortices
function total_velocity(vortices, strengths, point, self_index)
    total_v = [0.0, 0.0]
    for (j, vortex) in enumerate(vortices)
        if j == self_index    # skip self
            continue
        end
        r = compute_r(point, vortex)
        g = [0.0, 0.0, strengths[j]]
        v = induced_velocity(g, r)
        total_v += v[1:2]
    end
    return total_v
end

# Update vortex location
update_vortice(point, velocity, dt) = point + velocity * dt

# Leapfrog iteration step
function leapfrog(vortices, strengths, dt)
    new_vortices = copy(vortices)
    for (i, vortex) in enumerate(vortices)
        v = total_velocity(vortices, strengths, vortex, i)
        new_vortices[i] = update_vortice(vortex, v, dt)
    end
    return new_vortices
end

# Initialize vortices
vortices = [
    [0.0, -0.5],
    [0.0, 0.5],
    [1.0, -0.5],
    [1.0, 0.5],
]

# Circulation strengths (one per vortex!)
gamma = [1.0, -1.0, 1.0, -1.0]

dt = 0.01
steps = 4000

# Simulation loop
function simulate(vortices, strengths, dt, steps)
    xs = [zeros(steps) for _ in vortices]
    ys = [zeros(steps) for _ in vortices]

    for step in 1:steps
        for (i, vortex) in enumerate(vortices)
            xs[i][step] = vortex[1]
            ys[i][step] = vortex[2]
        end
        vortices = leapfrog(vortices, strengths, dt)
    end

    return xs, ys
end

xs, ys = simulate(vortices, gamma, dt, steps)

# Plot
plt = plot(xlabel="x", ylabel="y", title="Leapfrog Vortex Simulation", legend=:topright)
for i in 1:length(vortices)
    if i % 2 == 0
        color = :red
    else
        color = :blue
    end
    plot!(plt, xs[i], ys[i], label="Vortex $i", color=color)
    scatter!(plt, [xs[i][1]], [ys[i][1]], label="", marker=:circle, color=:black)  # starting points
end

display(plt)
println("Done")
