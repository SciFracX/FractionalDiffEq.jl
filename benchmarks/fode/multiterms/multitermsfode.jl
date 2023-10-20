using FractionalDiffEq, Plots, LinearAlgebra

h = 1e-2; tspan = (0, 30)
rightfun(x, y) = 172/125*cos(4/5*x)
multitermprob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [1, 4/5, -16/25, 0, 0, 0], tspan)

realfun(x)=sqrt(2)*sin(4*x/5+π/4)
sol1 = solve(multitermprob, h, PIEX())
sol2 = solve(multitermprob, h, PIIMRect())
sol3 = solve(multitermprob, h, PIIMTrap())
sol4 = solve(multitermprob, h, PIPECE())
#=
plot(sol1)
plot!(sol2)
plot!(sol3)
plot!(sol4)
=#
#plot!(collect(0:0.01:30), realfun.(collect(0:0.01:30)))

realsol = realfun.(collect(0:h:30))
err1 = norm(sol1.u-realsol)
err2 = norm(sol2.u-realsol)
err3 = norm(sol3.u-realsol)
err4 = norm(sol4.u-realsol)