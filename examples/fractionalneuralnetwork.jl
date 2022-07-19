using FractionalDiffEq, Plots

function sys!(du, u, p, t)
    du[1] = -0.05*u[2] - 0.05*u[3] + 0.01*tanh(u[2])
    du[2] = 0.05*u[1] + 0.02*u[2] + 0.01*tanh(u[1])
    du[3] = 0.1 - 0.2*u[3] + 0.05*u[1]*u[3] + 0.01*tanh(u[3])
end
prob = FractionalDifferenceSystem(sys!, 0.98, [1, -1, 0])
result = solve(prob, 997, GL())

plot(result[1, :], result[2, :], result[3, :], seriestype=:scatter)