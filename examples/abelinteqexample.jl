using FractionalDiffEq, Plots, SpecialFunctions
analytical(x)=exp(1+x)*erfc(sqrt(1+x))
tspan = LinRange(-1, 1, 100)
prob = FIEProblem([1, 1], [0.5, 0], 1, tspan)
sol = solve(prob, 20, SpectralUltraspherical())
solanalytical = analytical.(tspan)
plot(sol, title="Second kind Abel integral equation", label="Numerical")
plot!(tspan, solanalytical, ls=:dash, label="Analytical")