using FractionalDiffEq, Plots

a = 1;
mu = 4;

function Brusselator(du, u, p, t)
    du[1] = a - (mu + 1) * u[1] + u[1]^2 * u[2]
    du[2] = mu * u[1] - u[1]^2 * u[2]
    #du
end
alpha = [0.8; 0.8]
t0 = 0;
tfinal = 50;
y0 = [0.2; 0.03];
h = 0.01;
tspan = (t0, tfinal);
prob = FODESystem(Brusselator, alpha, y0, tspan)
#(t, y) = solve(prob, h, FLMMBDF())
#(t, y) = solve(prob, h, FLMMNewtonGregory())
(t, y) = solve(prob, h, FLMMTrap())

plot(y[1, :], y[2, :])
