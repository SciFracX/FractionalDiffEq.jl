"""
G2 Direct algorithm for fractional ordinary differential equations

```tex
@inproceedings{Guo2015FractionalPD,
  title={Fractional Partial Differential Equations and their Numerical Solutions},
  author={Boling Guo and Xueke Pu and Feng-Hui Huang},
  year={2015}
}
```
"""
struct G2Direct <: FractionalDiffEqAlgorithm end

"""
!!! info "Order of the problem"
    Please note the `L2Direct` method can only used for ``0 \\leq α \\leq 1``.
"""
struct L2Direct <: FractionalDiffEqAlgorithm end

function G2cₙⱼ(n, α, j)
    if j == 0
        return gamma(n-1-α)/(gamma(-α)*gamma(n))*(α^2/8 - α/4)
    elseif j == 1
        return gamma(n-2-α)/(gamma(-α)*gamma(n))*((1-α/2+α^2/8)*n+α^2/8-α/4-2)
    elseif 2 ≤ j ≤ n-1
        return gamma(n-j-1-α)/(gamma(-α)*gamma(n-j+2))*((n-j)^2 - (n-j)*α*(α+3)/2 + (α+1)*(α^3/8+α^2/2-1))
    elseif j == n
        return -α^3/8-α^2/2+1
    elseif j == n+1
        return (α^2+2*α)/8
    else
        return 0
    end
end

function solve(f, α, u0, T, h, ::G2Direct)
    n = Int64(floor(T/h))
    y = zeros(n+1)

    y[1]=u0
    for i in range(1, n, step=1)
        temp = zero(Float64)
        for j in range(1, n+1, step=1)
            temp += G2cₙⱼ(n, α, j-1)*y[i]
        end

        y[i+1] = h^α/G2cₙⱼ(n, α, i)*f(i*h, y[i]) - 1/G2cₙⱼ(n, α, i)*temp
    end
    
    return y
end

function solve(f, α, u0, T, h, ::L2Direct)
    n = Int64(floor(T/h))

    y = zeros(n+1)

    y[1]=u0
    for i in range(1, n-1, step=1)
        temp = zero(Float64)
        for j in range(1, i, step=1)
            temp += L2cₙⱼ(n, α, j-1)*y[i]
        end

        y[i+1] = h^α/L2cₙⱼ(n, α, i-1)*f(i*h, y[i]) - 1/L2cₙⱼ(n, α, i-1)*temp
    end
    
    return y
end
function L2cₙⱼ(n, α, j)
    temp = zero(Float64)
    if j == 0
        temp = bⱼ(α-1, n-1)
    elseif j == 1
        temp = -2*bⱼ(α-1, n-1) + bⱼ(α-1, n-2)
    elseif 2 ≤ j ≤ n-1
        temp = bⱼ(α-1, n-j+1)-2*bⱼ(α-1, n-j) + bⱼ(α-1, n-j-1)
    elseif j == n
        temp = 2^(2-α)-3
    elseif j == n+1
        temp = 1
    else
        temp = 0
    end
    return temp/gamma(3-α)
end
function bⱼ(α, j)
    return (j+1)^(1-α)-j^(1-α)
end
#=
h=0.5
T=5
fun(x, y) = 1 - y
result = solve(fun, 0.5, 0, T, h, G2Direct())
tspan=collect(0:h:T)

using Plots
plot(tspan, result)
=#