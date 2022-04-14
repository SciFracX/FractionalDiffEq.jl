using FractionalDiffEq, Plots
α=[0.94, 0.94, 0.94]; ϕ=[0.2, 0, 0.5]; τ=0.009; T=1.4; h=0.001
function delaychen(t, ϕ, y, k)
    a=35; b=3; c=27
    if k == 1
      return a*(y[2]-ϕ[1])
    elseif k == 2
      return (c-a)*ϕ[1]-y[1]*y[3]+c*y[2]
    elseif k == 3
      return y[1]*y[2]-b*ϕ[3]
    end
end
prob = FDDESystem(delaychen, ϕ, α, τ, T)
y, x=solve(prob, h, DelayABM())
plot(x[:, 1], x[:, 2], x[:, 3], title="Fractional Order Chen Delayed System")