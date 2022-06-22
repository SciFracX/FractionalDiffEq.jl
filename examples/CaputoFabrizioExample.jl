using FractionalDiffEq, Plots
a₁, a₂, a₃, a₄, a₅, a₆, a₇ = 3, 0.5, 4, 3, 4, 9, 4
function LotkaVolterra(du, u, p, t)
    du[1] = u[1]*(a₁-a₂*u[1]-u[2]-u[3])
    du[2] = u[2]*(1-a₃+a₄*u[1])
    du[3] = u[3]*(1-a₅+a₆*u[1]+a₇*u[2])    
end
u0 = [0.5, 0.9, 0.1]
tspan = (0, 50)
α = ones(3)*0.98;
prob = FODESystem(LotkaVolterra, α, u0, tspan)
sol = solve(prob, 0.01, NewtonPolynomial())

plot(sol[1, :], sol[2, :], sol[3, :])