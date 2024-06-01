using FractionalDiffEq, Plots

α = 1.8;
h = 1e-2;

# Analytical solution
analytical(x) = x .^ 1.8 .* mittleff(1.8, 2.8, -x .^ 1.8)
# Numerical solution
fun(x, y) = 1 - y
#prob = SingleTermFODEProblem(fun, 1.8, [0, 0], (0, 20))
prob = SingleTermFODEProblem(fun, α, [0, 0], (0, 20))
sol1 = solve(prob, h, PECE())
sol2 = solve(prob, h, GL())
sol3 = solve(prob, h, PIEX())

tspan = collect(0:0.01:20)
target = analytical.(tspan)

plot(sol1)
plot!(sol2)
plot!(sol3)
#plot!(tspan, target, lw=3, ls=:dash, label="Analytical")
