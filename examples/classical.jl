using FractionalDiffEq, Plots

function classical!(du, u, p, t)
    A, B = 1, 3
    du[1]=A-(B+1)*u[1]+u[1]^2*u[2]
    du[2]=B*u[1]-u[1]^2*u[2]
end
α = [0.8, 0.7]
u0 = [1.2, 2.8]
h=0.01
prob = FODESystem(classical!, α, u0)
result = solve(prob, h, 10, NonLinearAlg())
plot(result, vars=(1,2))