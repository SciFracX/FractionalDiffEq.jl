using FractionalDiffEq, Plots

a=1; mu=4
fdefun(t, y)=[a-(mu+1)*y[1]+y[1]^2*y[2]; mu*y[1]-y[1]^2*y[2]]
Jfdefun(t, y) = [-(mu+1)+2*y[1]*y[2] y[1]^2; mu-2*y[1]*y[2] -y[1]^2]
alpha=0.8
t0=0; tfinal=50; y0=[0.2; 0.03]
h=0.01
prob = FODESystem(fdefun, alpha, y0, t0, tfinal)
#(t, y) = solve(prob, Jfdefun, h, FLMMBDF())
(t, y) = solve(prob, Jfdefun, h, FLMMNewtonGregory())

plot(y[1, :], y[2, :])