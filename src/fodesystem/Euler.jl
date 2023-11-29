"""
The classical Euler method extended for fractional order differential equations.
"""
struct Euler <: FODESystemAlgorithm end

function solve(prob::FODEProblem, h, ::Euler)
    @unpack f, order, u0, tspan, p = prob
    t0=tspan[1]; tfinal=tspan[2]
    order = order[1]
    t = collect(Float64, t0:h:tfinal)
    N::Int = ceil(Int, (tfinal-t0)/h)
    l = length(u0)
    result = zeros(Float64, l, N+1)
    result[:, 1] = u0
    for n=1:N
        temp = zeros(Float64, l)
        tmp = zeros(Float64, l)
        for j=1:n
            f(tmp, result[:, j], p, t[j])
            temp = ((n-j+1)^order-(n-j)^order)*tmp
        end
        result[:, n+1] = result[:, 1] + h^order/gamma(order+1)*temp
    end
    return result
end