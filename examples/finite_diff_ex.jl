using FractionalDiffEq

K = 1
fdorder = 1.9
dx = pi/20
dt = 0.01
n = 2
xStart = 0
xEnd = pi

using Plots
plotlyjs()
U=solve(fdorder, dx, dt, xStart, xEnd, n, K, FiniteDiffEx())
plt=plot(x, t, U, st=:surface)