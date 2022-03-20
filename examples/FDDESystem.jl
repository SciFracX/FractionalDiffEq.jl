using FractionalDiffEq, Plots
α=[0.94, 0.94, 0.94]; ϕ=[0.2, 0, 0.5]; τ=0.009; T=1.4; h=0.001
function delaychen(x, y, z, xt, yt, zt, k)
    a=35; b=3; c=27
    if k == 1
      return a*(y-xt)
    elseif k == 2
      return (c-a)*xt-x*z+c*y
    elseif k == 3
      return x*y-b*zt
    end
end
prob = FDDESystem(delaychen, ϕ, α, τ, T)
x=solve(prob, h, DelayABMYuan())
plot(x[:, 1], x[:, 2], x[:, 3], title="Fractional Order Chen Delayed System")