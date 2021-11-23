using Plots

f(x)=exp(x)

T=30
h=0.05
tspan = collect(0.05:h:T)

result=f.(tspan)
plot(tspan, result)