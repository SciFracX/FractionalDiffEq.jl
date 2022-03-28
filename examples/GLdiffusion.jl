using FractionalDiffEq, Plots
plotlyjs()
rightfun(x, t) = -(1+x).*exp.(-t).*x.^3
T=1
M=100
N=M
h=1/M
tau=T/N
x=collect(0:h:1)
t=collect(0:tau:T)
alph=1.8
initial_condition(x) = x.^3
left_boundry(t) = 0
right_boundry(t) = exp.(-t)
d(x)  = gamma(2.2)*x.^2.8/6
(X, Y)=meshgrid(x, t)
sol=solve(alph, d, rightfun, M, N, initial_condition, left_boundry, right_boundry, GLDiff())
plot(X, Y, sol, st=:surface)