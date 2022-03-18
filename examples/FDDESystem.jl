using FractionalDiffEq
α=0.94; τ=0.9; T=1.4; h=0.1
function testf(x, y, z, xt, yt, zt, k)
    a=35
    b=3
    c=27
    if k == 1
      return a*(y-xt)
    elseif k == 2
      return (c-a)*xt-x*z+c*y
    elseif k == 3
      return x*y-b*zt
    end
end
(x, y, z)=solve(testf, α, τ, T, h, DelayABMYuan())

using Plots
plot(x, y)