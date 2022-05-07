using FractionalDiffEq, Plots

function testf!(t, y)
    p=[1 3]
    return [p[1]-(p[2]+1)*y[1] + y[1]^2*y[2]; p[2]*y[1]-y[1]^2*y[2]]
end

alpha = [0.8, 0.7]
t0=0; T=100
y0=[1.2; 2.8]
h=0.01
param=[1 3]
prob = FODESystem(testf!, alpha, y0, t0, T)
#(t, y) = solve(prob, h, PECE())
(t, y) = solve(prob, h, PIEX())
plot(y[1, :], y[2, :])
#=
plot(t, y[1, :])
plot!(t, y[2, :])
=#