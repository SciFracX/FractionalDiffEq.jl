using FractionalDiffEq, Plots
function chua!(du, x, p, t)
    a = 10.725; b = 10.593
    c = 0.268
    m0 = -1.1726
    m1 = -0.7872
    du[1] = a*(x[2]-x[1]-(m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))))
    du[2] = x[1]-x[2]+x[3]
    du[3] = -b*x[2]-c*x[3]
end
α = [0.93, 0.99, 0.92]
x0 = [0.2; -0.1; 0.1]
h = 0.01; tspan = (0, 10)
prob = FODESystem(chua!, α, x0, tspan)
result = solve(prob, h, NonLinearAlg())

plot(result[1, :], result[2, :], title="Chua System", legend=:bottomright)