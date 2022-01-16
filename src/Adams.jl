## Weird🤔

function testtestsolve(f, α, u0, h, T)
    N = Int64(floor(T/h))
    u = zeros(N+2)
    u[1] = u0
    #m = ceil(α)

    for j in range(0, N+1, step=1)
        u[j+1] = u0 + aⱼ(j, α, N)*f(j*h, u[j])*h^α/gamma(α+2)
    end

    return u
end

#=
function left(n, u0, h, m)
    temp = 0
    for j in range(0, m-1, step=1)
        temp += ((n+1)*h)^j/factorial(j)*u0^j
    end
    return temp
end
=#

function aⱼ(j, α, n)
    if j == 0
        return n^(α+1)-(n-α)*(n+1)^α
    elseif 1 ≤ j ≤ n
        return (n-j+2)^(α+1) - 2*(n-j+1)^(α+1) + (n-j)^(α+1)
    elseif j == n+1
        return 1
    end
end