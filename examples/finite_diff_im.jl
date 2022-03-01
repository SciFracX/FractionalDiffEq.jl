using FractionalDiffEq

K = 1
α = 0.5
dx = pi/2
dt = 0.5
n = 2
xStart = 0;
xEnd = pi    
x = collect(0:dx:xEnd)
t = collect(0:dt:n)
u0t = 0
uendt = 0
u0(x) = sin(x)

U=solve(α, dx, dt, xStart, xEnd, n, K, u0t, uendt, u0, FiniteDiffIm())

using Plots

plotlyjs()
plot(x, t, U, st=:surface)