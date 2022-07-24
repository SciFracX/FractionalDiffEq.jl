using FractionalDiffEq, Plots, BenchmarkTools

function test(du, u, p, t)
    du[1] = 1-4*u[1]+u[1]^2*u[2]
    du[2] = 3*u[1]-u[1]^2*u[2]
end

alpha = [0.8, 0.7]
t0=0; T=100
y0=[1.2; 2.8]
h=0.01
param=[1 3]
prob = FODESystem(test, alpha, y0, (t0, T))
sol = solve(prob, h, PIEX())
plot(sol, vars=(1, 2))