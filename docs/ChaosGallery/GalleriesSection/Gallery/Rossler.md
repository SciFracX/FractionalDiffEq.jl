using FractionalDiffEq

h=0.005
alpha = [0.9, 0.85, 0.95]
x0 = [0.5, 1.5, 0.1]
tf=120
function f(t, x, y, z, k)
    a, b, c = 0.5, 0.2, 10
    if k == 1
        return -y-z
    elseif k == 2
        return x+a*y
    elseif k == 3
        return b+z*(x-c)
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Rossler System")

![Rossler](./assets/Rossler.png)