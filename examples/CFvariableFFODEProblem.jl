using FractionalDiffEq, Plots

alpha=0.96;h=0.01;tfinal=100;
β(t) = 0.01+0.01*t
u0=[-0.2; 0.5; 0.2]
function fun(du, u, p, t)
    gama=10.814;lambda=14;a=0.2;b=0.15;
    du[1] = gama*(u[2]-a*sin(2*pi*b*u[1]))
    du[2] = u[1]-u[2]+u[3]
    du[3] = -lambda*u[2]
end

prob = FFODEProblem(fun, [alpha, β], u0, (1, tfinal))
result = solve(prob, h, AtanganaSeda())
plot3d(result[1, :], result[2, :], result[3, :])