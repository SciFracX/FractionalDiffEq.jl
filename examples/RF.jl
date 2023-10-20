using FractionalDiffEq, Plots

function RF(du, u, t, p)
    du[1] = u[2]*(u[3]-1+u[1]*u[1])+0.1*u[1];
    du[2] = u[1]*(3*u[3]+1-u[1]*u[1])+0.1*u[2];
    du[3] = -2*u[3]*(0.98+u[1]*u[2]);
end

prob = FODESystem(RF, [0.98, 0.98, 0.98], [0.1, 0.1, 0.1], (0, 100))
sol = solve(prob, 0.01, NonLinearAlg())
plot(sol, vars=(1,2,3))