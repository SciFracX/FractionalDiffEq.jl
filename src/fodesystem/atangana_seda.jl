function solve(prob::FODEProblem, ::AtanganaSedaAB; dt = 0.0)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack f, order, u0, tspan, p = prob
    order = order[1]
    t0 = tspan[1]; tfinal = tspan[2]
    t = collect(t0:dt:tfinal)
    N::Int = ceil(Int, (tfinal-t0)/dt)
    AB = 1-order+order/gamma(order)
    t = collect(Float64, t0:dt:tfinal)

    l = length(u0)
    # Initialization
    result = zeros(Float64, l, N+1)
    result[:, 1] = u0


    # Calculuate the first two value to start the computing
    tmp1=zeros(Float64, l); tmp2=zeros(Float64, l)
    f(tmp1, u0, p, t[1])
    result[:, 2] = u0+dt.*tmp1
    f(tmp2, result[:, 2], p, t[2])
    result[:, 3] = result[:, 2]+(dt/2).*(3*tmp2-tmp1)

    for n=3:N
        temp1 = zeros(l)
        temp2 = zeros(l)
        temp3 = zeros(l)

        temptemp1, temptemp2, temptemp3 = zeros(l), zeros(l), zeros(l)
        f(temp1, result[:, n-1], p, t[n-1])
        for j=3:n
            f(temp1, result[:, j-2], p, t[j-2])
            temptemp1 += ((n+1-j)^order-(n-j)^order)*temp1
            f(temp2, result[:, j-1], p, t[j-1])
            temptemp2 += (temp2-temp1)*((n+1-j)^order*(n-j+3+2*order)-(n-j)^order*(n-j+3+3*order))
            f(temp3, result[:, j], p, t[j])
            temptemp3 += (temp3-2*temp2+temp1)*((n-j+1)^order*(2*(n-j)^2+(3*order+10)*(n-j)+2*(order)^2+9*order+12)-(n-j)^order*(2*(n-j)^2+(5*order+10)*(n-j)+6*order^2+18*order+12))
        end
        temptemptemp = zeros(l)
        f(temptemptemp, result[:, n], p, t[n])

        result[:, n+1] = u0+(1-order)/AB*temptemptemp+((dt^order)*order/(AB*gamma(order+1)))*temptemp1 + ((dt^order)*order/(AB*gamma(order+2)))*temptemp2+((dt^order)*order/(2*AB*gamma(order+3)))*temptemp3
    end
    return FODESystemSolution(t, result)
end


"""
    solve(prob::FODESystem, dt, AtanganaSedaCF())

Atagana Seda method for Caputo-Fabrizio fractional order differential equations.

!!! tip
    Used for the Caputo Fabrizio fractional differential operators.

```tex
https://doi.org/10.1016/c2020-0-02711-8
```
"""
struct AtanganaSedaCF <: FODESystemAlgorithm end
#FIXME: Tests
function solve(prob::FODEProblem, ::AtanganaSedaCF; dt = 0.0)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    t=collect(Float64, t0:dt:tfinal)
    order=order[1]
    M=1-order+order/gamma(order)
    N=ceil(Int, (tfinal-t0)/dt)
    l=length(u0)
    result = zeros(l, N+1)
    temp1 = zeros(l); temp2 = zeros(l)
    f(temp1, u0, p, t[1])
    result[:, 2] = u0+temp1
    f(temp2, result[:, 2], p, t[2])
    result[:, 3] = result[:, 2] + (dt/2)*(3*temp2-temp1)
    tempn, tempn1, tempn2 = zeros(l), zeros(l), zeros(l)
    for n=3:N
        f(tempn, result[:, n], p, t[n])
        f(tempn1, result[:, n-1], p, t[n-1])

        f(tempn2, result[:, n-2], p, t[n-2])
        result[:, n+1] = result[:, n] + (1-order)/M*(tempn-tempn1)+order*M*dt*(23/12*tempn-4/3*tempn1+5/12*tempn2)
    end
    return FODESystemSolution(t, result)
end