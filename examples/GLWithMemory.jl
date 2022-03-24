using FractionalDiffEq
h=0.005; tf=30
alpha = [0.99, 0.99, 0.99]
x0 = [1, 0, 1]

function f!(du, u, p, t)
    a, b, c = 10, 28, 8/3
    du[1] = a*(u[2]-u[1])
    du[2] = u[1]*(b-u[3])-u[2]
    du[3] = u[1]*u[2]-c*u[3]
end
prob = FODESystem(f!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2])