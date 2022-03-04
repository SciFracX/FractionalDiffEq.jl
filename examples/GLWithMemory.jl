h=0.005
alpha = [0.99, 0.99, 0.99]
x0 = [1, 0, 1]
tf=30
function f(t, x, y, z, k)
    a, b, c = 10, 28, 8/3
    if k == 1
        return a*(y-x)
    elseif k == 2
        return x*(b-z)-y
    elseif k == 3
        return x*y-c*z
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2])