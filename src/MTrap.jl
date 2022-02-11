"""
Modified trapezoidal algorithms for system of fractional differential equations.

```tex
@Article{math8101675,
AUTHOR = {Zabidi, Nur Amirah and Abdul Majid, Zanariah and Kilicman, Adem and Rabiei, Faranak},
TITLE = {Numerical Solutions of Fractional Differential Equations by Using Fractional Explicit Adams Method},
DOI = {10.3390/math8101675}
}
```
"""
struct ModifiedTrap <: FractionalDiffEqAlgorithm end


function testsolve(f, α, x0, h, T, ::ModifiedTrap)
    equationsize = length(x0)
    N = Int64(floor(T/h))
    x = zeros(equationsize, N)

    # Initial conditions
    x[:, 1] = x0

    for i=1:equationsize
        for j=2:N
            x[i, j] = h^α[i]/gamma(α[i]+2)*((j-1)^(α[i]+1)-(j-α[i]-1)*j^α[i])*f(0, x0..., i) + x0[i] + h^α[i]/gamma(α[i]+2)*middle(i, j, h, α, f, x) + h^α[i]/gamma(α[i]+2)*f(h*j, genpara(i, j, h, α, f, x, equationsize)..., i)
        end
    end
    return x
end

function genpara(i, j, h, α, f, x, s)
    result=zeros(s)
    for k=1:s
        result[k]= x[k, j-1]+h^α[i]/gamma(α[i]+1)*f(h*(j-1), x[:, j-1]..., k)
    end
    return result
end

function middle(i, j, h, α, f, x)
    result=0
    for k=1:j-1
        result+=((j-k+1)^(α[i]+1)-2*(j-k)^(α[i]+1)+(j-k-1)^(α[i]+1))*f(k*h, x[:, k]..., i)
    end
    return result
end
#=
function f(t, x, y, z, k)
    if k == 1
        return 2*y^2
    elseif k == 2
        return t*x
    elseif k == 3
        return y*z
    end
end

function test(t, x, y, z, k)
    if k == 1
        return -x
    elseif k == 2
        return x-y^2
    elseif k == 3
        return y^2
    end
end

h=0.01
tspan = collect(h:h:1)
result = testsolve(test, [0.99, 0.99, 0.99], [1, 0, 0], h, 1)
using Plots

plot(tspan, result[2, :])
=#
#=
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
h=0.01
T=50
x0=[0.1, 0.2, 0.3]
#prob=FODESystem(qi, alpha, x0)
result = testsolve(qi, alpha, x0, h, T)
using Plots
plot(result[1, :], result[2, :])
=#