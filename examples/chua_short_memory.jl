using FractionalDiffEq, Plots

function chua(du, u, p, t)
    a=10.725; b=10.593; c=0.268
    m0=-1.1726; m1=-0.7872
    du[1] = a*(u[2]-u[1]-m1*u[1]+0.5*(m0-m1)*(abs(u[1]+1)-abs(u[1]-1)))
    du[2] = u[1]-u[2]+u[3]
    du[3] = -b*u[2]-c*u[3]
end

α = [0.93, 0.99, 0.92];
u0 = [0.2; -0.1; 0.1];
h = 0.001;
tspan = (0, 500);
prob = FODESystem(chua, α, u0, tspan)
memory = 10000
result = solve(prob, h, tn, NonLinearAlg(), memory)
plot(result[:, 1], result[:, 2], title="Chua System with Short Memory Effect", legend=:bottomright)