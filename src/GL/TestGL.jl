import FractionalDiffEq: FractionalDiffEqAlgorithm, SingleTermFODEProblem, solve

struct GL <: FractionalDiffEqAlgorithm end

# Some points are a little big
function solve(prob::SingleTermFODEProblem, u0, T, ::GL)
    f, α, h = prob.f, prob.α, prob.h
    N = Int64(floor(T/h))

    c = zeros(N)

    cp = 1
    for j = 1:N
        c[j] = (1-(1+α)/j)*cp
        cp = c[j]
    end

    # Initialization
    y = zeros(N)
    y[1] = u0

    for i = 2:N
        y[i] = f(y[i-1])*h^α - Cq(y, c, i)
    end

    return y
end

function Cq(r, c, k)
    temp = 1
    for j = 1:k-1
        temp = temp + c[j]*r[k-j]
    end
    return temp
end