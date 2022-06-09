"""
    solve(prob::SingleTermFODEProblem, h, AS())

Atangana-Seda method for Caputo single term FODE.
"""
struct AtanganaSeda <: FractionalDiffEqAlgorithm end

function solve(prob::SingleTermFODEProblem, h, ::AtanganaSeda)
    @unpack f, α, u0, tspan = prob
    t0 = tspan[1]; tfinal = tspan[2]
    N = ceil(Int, (tfinal-t0)/h)
    t = collect(t0:h:tfinal)
    u=zeros(N+1)

    u[1]=u0
    u[2]=u[1]+h*f(t[1], u[1]) # One step Euler method to first evaluate the second value

    u[3]=u[2]+(h/2)*(3*f(t[2], u[2])-f(t[1], u[1])) # Two-step Adams-Bashforth method to evaluate the third value

    for n=3:N
        temp1=0
        temp2=0
        temp3=0
        for j=3:n
            temp1 += f(t[j-2], u[j-2])*((n-j+1)^α-(n-j)^α)
            temp2 += (f(t[j-1], u[j-1])-f(t[j-2], u[j-2]))*((n-j+1)^α*(n.-j.+3 .+2*α).-(n.-j).^α.*(n.-j.+3 .+3*α))
            temp3 += ((f(t[j], u[j])-2*f(t[j-1], u[j-1])+f(t[j-2], u[j-2])))*(((n-j+1)^α*(2*(n-j)^2+(3*α+10)*(n-j)+2*α^2+9*α+12)) - (n-j)^α*(2*(n-j)^2+(5*α+10)*(n-j)+6*α^2 +18*α+12))
        end
        u[n+1] = @. u[1]+(h^α/gamma(α+1))*temp1+(h^(α)/gamma(α+2))*temp2 + (h^α/(2*gamma(α+3)))*temp3
    end
    return FODESolution(t, u)
end