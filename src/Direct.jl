import FractionalDiffEq.FractionalDiffEqAlgorithm

struct G2Direct <: FractionalDiffEqAlgorithm end

function G2cₙⱼ(n, α, j)
    if j == 0
        return gamma(n-1-α)/(gamma(-α)*gamma(n))*(α^2/8-α/4)
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


    y = zeros(n)
    temp = zero(Float64)
    y[1]=u0
    for i in range(1, n-1, step=1)
        


        for j in range(1, n-1, step=1)
            temp += G2cₙⱼ(n, α, j)*y[i]
        end

        y[i+1] = h^α/G2cₙⱼ(n+1, α, n)*f(i*h, y[i]) - 1/G2cₙⱼ(n+1, α, n)*temp
    end
    
    return y
end


function L2cₙⱼ(n, α, j)
    if j == 0
        return 1
    end
end

function bⱼ(n, α, j)
    return (j+1)^(1-α)-j^(1-α)
end
