using FractionalDiffEq, Plots

# Analytical solution
analytical(x) = x .^ 1.8 .* mittleff(1.8, 2.8, -x .^ 1.8)
# Numerical solution
fun(y, p, t) = 1 - y
prob = SingleTermFODEProblem(fun, 1.8, [0, 0], (0, 20))
sol = solve(prob, 0.01, PECE())
tspan = collect(0:0.01:20)
target = analytical.(tspan)

plot(sol, linewidth = 5, label = "Numerical", legend = :bottomright)
plot!(tspan, target, lw = 3, ls = :dash, label = "Analytical")
