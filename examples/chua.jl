using FractionalDiffEq, Plots
function chua!(du, x, p, t)
    a, b, c, m0, m1 = p
    du[1] = a *
            (x[2] - x[1] - (m1 * x[1] + 0.5 * (m0 - m1) * (abs(x[1] + 1) - abs(x[1] - 1))))
    du[2] = x[1] - x[2] + x[3]
    du[3] = -b * x[2] - c * x[3]
end
α = [0.93, 0.99, 0.92]
x0 = [0.2; -0.1; 0.1]
h = 0.001;
tspan = (0, 50);
p = [10.725, 10.593, 0.268, -1.1726, -0.7872]
prob = FODESystem(chua!, α, x0, tspan, p)
sol = solve(prob, h, NonLinearAlg())

plot(sol, vars = (1, 2), title = "Fractional Chua System", legend = :bottomright)
