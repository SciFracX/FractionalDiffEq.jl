using FractionalDiffEq
using Plots

tmp = solve(0.7, 1.8, 1, 21, 148, FPDEMatrixDiscrete())



YS = reshape(tmp, 19, 147)
YS = reverse(YS, dims=2)
U = YS


rows, columns = size(U)

U = [zeros(1, columns); U; zeros(1, columns)]
U=[zeros(1, 21)' U]

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)

plotly()

plot(XX, YY, U, st=:surface)