using FractionalDiffEq, Plots

e(x)=1+0*x
f(x)=0*x
xx = LinRange(-1, 1, 100)
sol = solve(f, e, 20, xx, SpectralUltraspherical())
plot(xx, sol)