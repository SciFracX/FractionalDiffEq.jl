using SpecialFunctions

"""
@article{
title={A predictor-corrector approach for the numerical solution of fractional differential equations},
author={Diethelm, Kai and Ford, Neville J. and Freed, Alan D.}
doi={https://doi.org/10.1023/A:1016592219341}
}
"""

function solve(f, α, u0, T, h)
    N=T/h
    y=zeros(Int64(N+1))

    for n in range(0, N, step=1)
        y[Int64(n+1)]=u0 + h^α/gamma(α+2)*f((n+1)*h, predictor(f, y, α, n, h, u0))+h^α/gamma(α+2)+right(f, y, α, n, h)
    end

    return y
end

function right(f, y, α, n, h)
    temp = 0

    for j in range(0, n, step=1)
        temp+=A(j, n, α)*f(j*h, y[Int64(j+1)])
    end

    return h^α/gamma(α+2)*temp
end

function predictor(f, y, α, n, h, u0)
     predict = 0

     for j in range(0, n, step=1)
        predict+=B(j, n, α, h)*f(j*h, y[Int64(j+1)])
     end

     return u0+predict
end


function A(j, n, α)
    if j==0
        return n^(α+1)-(n-α)*(n+1)^α
    elseif 1 ≤ j ≤ n
        return (n-j+2)^(α+1)+(n-j)^(α+1)-2*(n-j+1)^(α+1)
    elseif j==n+1
        return 1
    end
end

function B(j, n, α, h)
    return h^α/α*((n+1-j)^α-(n-j)^α)
end

fun(x, y) = 1-y
result=solve(fun, 0.75, 0, 1, 0.01)
tspan=collect(0:0.01:1)
print(result)
#plot(tspan, result)