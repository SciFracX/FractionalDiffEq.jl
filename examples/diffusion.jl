using FractionalDiffEq
using Plots

α = 0.7
β = 1.8
κ = 1
T = 2
m = 21
n = 148

result = solve(α, β, κ, T, m, n, FPDEMatrixDiscrete())

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)

plotlyjs()

plot(XX, YY, result, st=:surface)