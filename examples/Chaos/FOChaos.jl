#=======Fractional Order Liu System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [0.2, 0, 0.5]
tf=100
function f(t, x, y, z, k)
    a, b, c, e, n, m = 1, 2.5, 5, 1, 4, 4
    if k == 1
        return -a*x-e*y^2
    elseif k == 2
        return b*y-n*x*z
    elseif k == 3
        return -c*z+m*x*y
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Liu System")


#=======Fractional Order Duffing System=======#
using FractionalDiffEq

h=0.005
alpha = [0.9, 1]
x0 = [0.21, 0.31]
tf=100
function f(t, x, y, k)
    α, δ, ω = 0.15, 0.3, 1
    if k == 1
        return y
    elseif k == 2
        return x-α*y-x^3+δ*cos(ω*t)
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Duffing System")


#=======Fractional Order Chen System=======#
using FractionalDiffEq

h=0.005
alpha = [0.9, 0.9, 0.9]
x0 = [-9, -5, 14]
tf=100
function f(t, x, y, z, k)
    a, b, c, d = 35, 3, 28, -7
    if k == 1
        return a*(y-x)
    elseif k == 2
        return d*x-x*z+c*y
    elseif k == 3
        return x*y-b*z
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Chen System")


#=======Fractional Order Van der Pol Oscillator=======#
using FractionalDiffEq

h=0.005
alpha = [1.2, 0.8]
x0 = [0.2, -0.2]
tf=60
function f(t, x, y, k)
    ϵ = 1
    if k == 1
        return y
    elseif k == 2
        return -x-ϵ*(x^2-1)*y
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Van der Pol Oscillator")


#=======Fractional Order Lotka Volterra System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [1, 1.4, 1]
tf=50
function f(t, x, y, z, k)
    a, b, c, d, e, p, s = 1, 1, 1, 1, 2, 3, 2.7
    if k == 1
        return a*x+e*x^2-b*x*y-s*z*x^2
    elseif k == 2
        return -c*y+d*x*y
    elseif k == 3
        return -p*z+s*z*x^2
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lotka Volterra System")


#=======Fractional Order Lu System=======#
using FractionalDiffEq

h=0.005
alpha = [0.985, 0.99, 0.98]
x0 = [0.2, 0.5, 0.3]
tf=60
function f(t, x, y, z, k)
    a, b, c = 36, 3, 20
    if k == 1
        return a*(y-x)
    elseif k == 2
        return -x*z+c*y
    elseif k == 3
        return x*y-b*z
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lu System")


#=======Fractional Order Rossler System=======#
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


#=======Fractional Order Arneodo System=======#
using FractionalDiffEq

h=0.005
alpha = [0.97, 0.97, 0.96]
x0 = [-0.2, 0.5, 0.2]
tf=100
function f(t, x, y, z, k)
    b1, b2, b3, b4 = -5.5, 3.5, 0.8, -1.0
    if k == 1
        return y
    elseif k == 2
        return z
    elseif k == 3
        return -b1*x-b2*y-b3*z+b4*x^3
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Arneodo System")


#=======Fractional Order Genesio-Tesi System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [-0.1, 0.5, 0.2]
tf=100
function f(t, x, y, z, k)
    b1, b2, b3, b4 = 1, 1.1, 0.4, 1.0
    if k == 1
        return y
    elseif k == 2
        return z
    elseif k == 3
        return -b1*x-b2*y-b3*z+b4*x^2
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Genesio-Tesi System")