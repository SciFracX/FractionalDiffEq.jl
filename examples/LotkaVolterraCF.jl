using FractionalDiffEq

function lt(du, u, p, t)
    a1, a2, a3, a4, a5, a6, a7 = 3, 0.5, 4, 3, 4, 9, 4
    du[1] = u[1] * (a1 - a2 * u[1] - u[2] - u[3])
    du[2] = u[2] * (1 - a3 + a4 * u[1])
    du[3] = u[3] * (1 - a5 + a6 * u[1] + a7 * u[2])
end

α = [0.98; 0.98; 0.98];
u0 = [0.5; 0.9; 0.1]
prob = FODESystem(lt, α, u0, (0, 200))
sol = solve(prob, 0.01, NewtonPolynomial())

using Plots

plot3d(sol[1, :], sol[2, :], sol[3, :])
