using FractionalDiffEq
using Plots

result = solve(0.7, 1.8, 2, 21, 148, FPDEMatrixDiscrete())

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)

plotlyjs()

plot(XX, YY, result, st=:surface)
savefig("diffusion.html")