function solve(prob::FODESystem, h, ::Euler)
    @unpack f, α, u0, tspan, p = prob
    t0=tspan[1]; tfinal=tspan[2]
    α = α[1]
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
            temp = ((n-j+1)^α-(n-j)^α)*tmp
        end
        result[:, n+1] = result[:, 1] + h^α/gamma(α+1)*temp
    end
    return result
end