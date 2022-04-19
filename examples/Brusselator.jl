using FractionalDiffEq, Plots

α = [0.8, 0.8]
u0 = [0.2, 0.03]
h = 0.001

function Brusselator!(du, u, p, t)
    a, μ = 1, 4
    du[1] = a-(μ+1)*u[1]+(u[1])^2*u[2]
    du[2] = μ*u[1]-(u[1])^2*u[1]
end

prob = FODESystem(Brusselator!, α, u0, 100)
result = solve(prob, h, GL())

# Phase plane
plot(result[:, 1], result[:, 2])

# Time plane
plot(collect(0:h:100), result[:, 1])
plot!(collect(0:h:100), result[:, 2])