# FIXME: Not well done!

"""
http://www.koreascience.or.kr/article/JAKO200833338752380.pdf
"""

import FractionalDiffEq.solve

using Plots, SpecialFunctions

using MittagLeffler

# Fractional Euler method

"""
Ahmed, H. (2018). FRACTIONAL EULER METHOD; AN EFFECTIVE TOOL FOR SOLVING FRACTIONAL DIFFERENTIAL EQUATIONS. Journal of the Egyptian Mathematical Society, 26(1), 38-43. doi: 10.21608/JOEMS.2018.9460
"""

function solve(f, α, u0, h, T)
    N = Int64(floor(T/h))
    #m=Int64(ceil(α))
    
    y=ones(N+1)
    y[1]=u0

    for i in range(1, N, step=1)
        #y[i+1]=left(i, u0, h, m) + h^α/gamma(α+1)*coeff(i, α, N)*f(i*h, y[i])
        y[i+1] = y[i] + h^α/gamma(α+1)*f(i*h, y[i])
    end

    return y
end

function coeff(j, α, n)
    return (n-j+1)^α-(n-j)^α
end

function left(n, u0, h, m)
    temp = 0
    for j in range(0, m, step=1)
        temp += ((n+1)*h)^j/factorial(i)*u0^i
    end
    return temp
end

fun(t, y) = y

result=solve(fun, 0.5, 0, 0.001, 5)
tspan=collect(0:0.001:5)

plot(tspan, result)

#Analytical solution
target = []


#MittagLeffler.jl doesn't support array argument
for i in 0:0.001:5
    push!(target, mittleff(0.5, -i^0.5))
end

plot!(tspan, target)