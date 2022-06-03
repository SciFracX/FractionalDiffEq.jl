function solve(FDDE::FDDEProblem, h, ::DelayPECE)
    @unpack f, ϕ, α, τ, tspan = FDDE
    t = collect(0:h:tspan)
    maxn = length(t)
    yp = zeros(maxn)
    y = copy(t)
    y[1] = ϕ(0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+b(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp)...)
        end
        yp[n+1] = yp[n+1]/gamma(α)+ϕ(0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+a(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp)...)
        end

        y[n+1] = y[n+1]/gamma(α)+h^α*f(t[n+1], yp[n+1], v(ϕ, n+1, τ, h, y, yp)...)/gamma(α+2) + ϕ(0)
    end

    V = []
    for n = 1:maxn
        push!(V, v(ϕ, n, τ, h, y, yp))
    end

    delayed = zeros(length(τ), length(V))
    for i=1:length(V)
        delayed[:, i] = V[i]
    end

    
    return delayed, y
end

function a(j, n, α, h)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(α+1)-(n-α)*(n+1)^α
    elseif j == n
        result = 2*(2^(α+1)-1)
    else
        result = (n-j+2)^(α+1) + (n-j)^(α+1) - 2*(n-j+1)^(α+1)
    end
    return result*h^α / (α*(α + 1))
end

function b(j, n, α, h)
    return h^α/α*((n-j+1)^α - (n-j)^α)
end

function v(ϕ, n, τ, h, y, yp)
    if maximum(τ) > n*h
        return ϕ.((n-1)*h.-τ)
    else
        m = floor.(Int, τ./h)
        δ = m.-τ./h

        function judge(m)
            temp1 = findall(x->x>1, m)
            temp2 = findall(x->x==1, m)#FIXME: Another case for x == 1

            result = zeros(length(m))
            if length(temp1) == length(m)
                for i=1:length(m)
                    result[i] = δ[i]*y[n-m[i]+2]+(1-δ[i])*y[n-m[i]+1]
                end
                return result
            end
        end
        return judge(m)
    end
end