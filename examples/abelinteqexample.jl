using FractionalDiffEq, Plots, SpecialFunctions
analytical(x)=exp(1+x)*erfc(sqrt(1+x))
xx = LinRange(-1, 1, 100)
prob = FIEProblem([1, 1], [1, 0.5], 1)
sol = solve(prob, 20, xx, SpectralUltraspherical())
solanalytical = analytical.(xx)
plot(xx, sol, title="Second kind Abel integral equation", label="Numerical")
plot!(xx, solanalytical, ls=:dash, label="Analytical")