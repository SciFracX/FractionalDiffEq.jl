using FractionalDiffEq, Plots, SpecialFunctions
analytical(x)=exp(1+x)*erfc(sqrt(1+x))
e(x)=1+0*x
f(x)=0*x
xx = LinRange(-1, 1, 100)
sol = solve(f, e, 20, xx, SpectralUltraspherical())
solanalytical = analytical.(xx)
plot(xx, sol, title="Second kind Abel integral equation", label="Numerical")
plot!(xx, solanalytical, ls=:dash, label="Analytical")