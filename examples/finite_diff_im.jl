using FractionalDiffEq

K = 1
α = 0.5
dx = pi/20
dt = 0.005
n = 2
xStart = 0;
xEnd = pi    
x = collect(0:dx:xEnd)
t = collect(0:dt:n)

U=solve(fdorder, dx, dt, xStart, xEnd, n, K, FiniteDiffEx())

using Plots

plotlyjs()
plot(x, t, U, st=:surface)