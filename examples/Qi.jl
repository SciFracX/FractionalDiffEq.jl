using FractionalDiffEq
using Plots

function qi(t, x, y, z, k)
    a, b, c, d, r = 35, 8/3, 80, -1, 1
    if k == 1
        return -a*x+a*y+r*y*z
    elseif k == 2
        return c*x+d*y-x*z
    elseif k == 3
        return -b*z+x*y
    end
end

alpha = [0.98, 0.98, 0.98]
h=0.001
T=50
x0=[0.1, 0.2, 0.3]
prob=FODESystem(qi, alpha, x0)
x, y, z = solve(prob, h, T, GLWithMemory())

plot(x, y)