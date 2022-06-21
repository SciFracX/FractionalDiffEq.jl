using FractionalDiffEq, Plots
α = [0.9; 0.9; 0.9];
alpha = [0.9, 0.8, 1]
u0 = [1-0.001; 0.001; 0];
tspan = (0, 100);
h = 0.01;
function SIR!(du, u, p, t)
    β, γ = 0.4, 0.04
    du[1] = -β*u[1]*u[2]
    du[2] = β*u[1]*u[2]-γ*u[2]
    du[3] = γ*u[2]
end

prob = FODESystem(SIR!, alpha, u0, tspan)
t, sol1 = solve(prob, h, FLMMBDF())
t, sol2 = solve(prob, h, FLMMNewtonGregory())
t, sol3 = solve(prob, h, FLMMTrap())
t, sol4 = solve(prob, h, PECE())
t, sol5 = solve(prob, h, PIEX())
sol6 = solve(prob, h, NonLinearAlg())
sol7 = solve(prob, h, GL())

plot(t, sol2[1, :])
plot!(t, sol2[2, :])
plot!(t, sol2[3, :])