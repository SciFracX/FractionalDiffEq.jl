using FractionalDiffEq, SpecialFunctions, Plots
tfinal = 50;
h = 0.01;
α = 0.98 * ones(3);
tspan = (0, 200)
a = 1;
b = 0.1;
function fun(du, u, p, t)
    du[1] = -sin(u[2]) - b * u[1]
    du[2] = -sin(u[3]) - b * u[2]
    du[3] = -sin(u[1]) - b * u[3]
end

u0 = [-1, 1, 1]

prob = FODESystem(fun, α, u0, tspan)
sol = solve(prob, h, AtanganaSedaCF())
plot(sol[1, :], sol[2, :], sol[3, :])
