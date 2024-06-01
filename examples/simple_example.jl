using FractionalDiffEq, Plots
fun(x, y) = 1 - y
prob = SingleTermFODEProblem(fun, 0.5, 0, (0, 5))
sol = solve(prob, 0.01, PECE())
plot(sol, linewidth = 2, legend = :bottomright)
