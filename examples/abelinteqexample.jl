using FractionalDiffEq

e(x)=1+0*x
f(x)=0*x
xx = LinRange(-1, 1, 100)
result=solve(f, e, 20)
sol=myeval(result, xx, 2)

using Plots
plot(xx, sol)