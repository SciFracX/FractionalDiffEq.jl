using FractionalDiffEq, Plots
α=1;β=1;h=0.01;tfinal=50
u0 = [-2, 1, -1]
a=10;b=28;c=8/3
function fun(du, u, p, t)
    a=10;b=28;c=8/3
    du[1] = a*(u[2]-u[1])
    du[2] = (b-u[3])*u[1]-u[2]
    du[3] = u[1]*u[2]-c*u[3]
end
prob = FFODESystem(fun, [α, β], u0, (0, tfinal))
sol = solve(prob, h, AtanganaSeda())
plot3d(sol[1, :], sol[2, :], sol[3, :], title="Fractal-fractional Order Lorenz System")