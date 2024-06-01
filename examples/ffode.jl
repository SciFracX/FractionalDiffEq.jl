using FractionalDiffEq, Plots

fun(t, y) = -3 * t + 17
u0 = -0.1
h = 0.01;
α = 0.75;
β = 0.5;
tfinal = 1;
tspan = (0, tfinal);
prob = FFMODEProblem(fun, [α, β], u0, tspan)
sol = solve(prob, h, AtanganaSeda())

plot(sol)
