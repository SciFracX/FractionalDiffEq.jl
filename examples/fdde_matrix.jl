using FractionalDiffEq

limit=100
t0=0
T=1
tau=3.1416
h=0.5
alpha=0.4
function x0(t)
    return [sin(t)*cos(t); sin(t)*cos(t); cos(t)^2-sin(t)^2; cos(t)^2-sin(t)^2]
end
A=[0 0 1 0; 0 0 0 1; 0 -2 0 0; -2 0 0 0]
B=[0 0 0 0; 0 0 0 0 ;-2 0 0 0; 0 -2 0 0]
f=[0; 0; 0; 0]

result=solve(limit, alpha, A, B, f, t0, x0, T, tau, h, MatrixForm())

using Plots

plot(result[:, 1], result[:, 3])