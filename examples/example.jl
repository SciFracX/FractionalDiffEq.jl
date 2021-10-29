using FractionalDiffEq
using Plots
using MittagLeffler

#Analytical solution
target = []

for i in 0:0.01:20
    push!(target, i^1.8*mittleff(1.8,2.8,-i^1.8))
end

#Numerical solution
fun(x, y) = 1-y
prob = FDEProblem(fun, 1.8, 0, 20, 0.01)
result=solve(prob)
tspan=collect(0:0.01:20)

#Want to use a more elegant plotting backend, especially for the LaTeX rendering, failed anyway (T_T)
#pgfplotsx()
gr()

plot(tspan, result, title="D^1.8 y(x)=1-y, y(0)=0", linewidth=5, label="Numerical", legend=:bottomright)
plot!(tspan, target, lw=3, ls=:dash, label="Analytical")