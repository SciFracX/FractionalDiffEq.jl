using FractionalDiffEq, Plots
gr()
h=0.1; alpha=[1.8, 1.8, 1.8]; phi=[0.01, 0.02]; tau=0.2
function test(du, u, ϕ, t)
    du[1] = -0.8*u[1]  -0.5*ϕ[1] - u[1]*ϕ[2]
    du[2] = -0.5*u[2] + 0.3*ϕ[2] - u[2]*ϕ[1]
end
prob = FDDESystem(test, phi, alpha, tau, 100)
sol = solve(prob, h, DelayABM())
plot(sol, vars=(1, 2))