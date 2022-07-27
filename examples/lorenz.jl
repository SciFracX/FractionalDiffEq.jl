using FractionalDiffEq, Plots
function lorenz(du, u, p, t)
    a, b, c, d = p
    du[1] = a*(u[2]-u[1])
    du[2] = c*u[1]-u[1]*u[3]+d*u[2]
    du[3] = u[1]*u[2]-b*u[3]
end
α0 = [0.96, 0.96, 0.96]
x0 = [1, 2, 3]; tspan=(0, 20)
h=0.01; p=[40, 3, 10, 15]
prob=FODESystem(lorenz, α0, x0, tspan, p)

sol=solve(prob, h, GL())
plot(sol, vars=(1, 3))