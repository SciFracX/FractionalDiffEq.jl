struct GL <: FractionalDiffEqAlgorithm end

# Need to be verified
function solve(prob::SingleTermFODEProblem, T, b, ::GL)
    f, α, h = prob.f, prob.α, prob.h
    n = ceil(Int, T/h)
    y = zeros(n)# Prelocate result
    y[1] = 0

    for k in range(2, n, step=1)
        y[k]=-b*h^α*y[k-1]-middle(k, y, α)+h^α*f(h*k)
    end

    return y
end

function middle(k, y, α)
    temp = 0

    for j in range(1, k, step=1)
        temp += coeff(α, j)*y[k-j+1]
    end
    return temp
end
#This can be computed using FFTW
function coeff(α, k)
    if k==0
        return 1
    else
        return (1-(α+1)/k)*coeff(α, k-1)
    end
end