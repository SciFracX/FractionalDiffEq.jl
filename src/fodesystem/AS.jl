"""
    solve(prob::FODESystem, h, AtanganaSedaAB())

Solve Atangana-Baleanu fractional order differential equations using Newton Polynomials.
"""
struct AtanganaSedaAB <: AbstractFDEAlgorithm end

function solve(prob::FODESystem, h, ::AtanganaSedaAB)
    @unpack f, α, u0, tspan = prob
    α = α[1]
    t0 = tspan[1]; tfinal = tspan[2]
    N = ceil(Int, (tfinal-t0)/h)
    AB = 1-α+α/gamma(α)
    t = collect(t0:h:tfinal)

    l = length(u0)
    # Initialization
    result = zeros(l, N+1)
    result[:, 1] = u0


    # Calculuate the first two value to start the computing
    tmp1=zeros(l); tmp2=zeros(l)
    f(tmp1, u0, nothing, t[1])
    result[:, 2] = u0+h.*tmp1
    f(tmp2, result[:, 2], nothing, t[2])
    result[:, 3] = result[:, 2]+(h/2).*(3*tmp2-tmp1)

    for n=3:N
        temp1 = zeros(l)
        temp2 = zeros(l)
        temp3 = zeros(l)

        temptemp1, temptemp2, temptemp3 = zeros(l), zeros(l), zeros(l)
        f(temp1, result[:, n-1], nothing, t[n-1])
        for j=3:n
            f(temp1, result[:, j-2], nothing, t[j-2])
            temptemp1 += ((n+1-j)^α-(n-j)^α)*temp1
            f(temp2, result[:, j-1], nothing, t[j-1])
            temptemp2 += (temp2-temp1)*((n+1-j)^α*(n-j+3+2*α)-(n-j)^α*(n-j+3+3*α))
            f(temp3, result[:, j], nothing, t[j])
            temptemp3 += (temp3-2*temp2+temp1)*((n-j+1)^α*(2*(n-j)^2+(3*α+10)*(n-j)+2*(α)^2+9*α+12)-(n-j)^α*(2*(n-j)^2+(5*α+10)*(n-j)+6*α^2+18*α+12))
        end
        temptemptemp = zeros(l)
        f(temptemptemp, result[:, n], nothing, t[n])

        result[:, n+1] = u0+(1-α)/AB*temptemptemp+((h^α)*α/(AB*gamma(α+1)))*temptemp1 + ((h^α)*α/(AB*gamma(α+2)))*temptemp2+((h^α)*α/(2*AB*gamma(α+3)))*temptemp3
    end
    return result
end