using FractionalDiffEq

function lorenz(t, x, y, z, k)
    a=40
    b=3
    c=10
    d=15
    if k==1
        return a*(y-x)
    elseif k==2
        return c*x-x*z+d*y
    elseif k==3
        return x*y-b*z
    end
end

α0 = [0.96, 0.96, 0.96]
x0 = [1, 2, 3]
h=0.001
prob=FODESystem(lorenz, α0, x0)
T=20
x, y, z=solve(prob, h, T, GLWithMemory())
using Plots
plot(x, z)

####################################
# Or use the detailed model in FractionalDiffEq.jl

a=40; b=3; c=10; d=15
prob = FractionalLorenz(a, b, c, d, 0.96)
x, y, z = solve(prob, h, T, LorenzADM())

using Plots
plot(x, z)