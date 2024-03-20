@concrete mutable struct ABMCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    constant_algs
    p

    x
    x0
    x1
    N
    Ndelay

    dt
    kwargs
end

Base.eltype(::ABMCache{iip, T}) where {iip, T} = T

"""
    solve(FDDE::FDDEProblem, dt, DelayABM())

Use the Adams-Bashforth-Moulton method to solve fractional delayed differential equations.

### References

```tex
@inproceedings{Bhalekar2011APS,
  title={A PREDICTOR-CORRECTOR SCHEME FOR SOLVING NONLINEAR DELAY DIFFERENTIAL EQUATIONS OF FRACTIONAL ORDER},
  author={Sachin Bhalekar and Varsha Daftardar-Gejji},
  year={2011}
}
```
"""
struct DelayABM <: FDDEAlgorithm end
#FIXME: There are still some improvments about initial condition
#FIXME: Fix DelayABM method for FDDESystem : https://www.researchgate.net/publication/245538900_A_Predictor-Corrector_Scheme_For_Solving_Nonlinear_Delay_Differential_Equations_Of_Fractional_Order
#TODO: Need more works
function SciMLBase.__init(prob::FDDEProblem, alg::DelayABM; dt=0.0, kwargs...)
    @unpack f, order, u0, h, tspan, p, constant_lags = prob
    τ = constant_lags[1]
    T = eltype(u0)
    l = length(u0)
    iip = SciMLBase.isinplace(prob)
    t0 = tspan[1]; tfinal = tspan[2]
    N = round(Int, (tfinal-t0)/dt)
    mesh = collect(T, t0:dt:tfinal)
    Ndelay = round(Int, τ/dt)
    #x1 = zeros(Float64, Ndelay+N+1)
    #x = zeros(Float64, Ndelay+N+1)
    #x1 = [zeros(T, l) for i=1:(Ndelay+N+1)]
    #x = [zeros(T, l) for i=1:(Ndelay+N+1)]
    x1 = _generate_similar_array(u0, Ndelay+N+1, u0)
    x = _generate_similar_array(u0, Ndelay+N+1, u0)
    #x1[Ndelay+N+1] = 0
    
    #x[Ndelay+N+1] = 0

    # History function handling
    #FIXME: When the value of history function h is different with the initial value?
    println(x[1])
    x[Ndelay] = length(u0) == 1 ? ([h(p, 0)]) : (h(p, 0))
    for i=1:Ndelay
        x[i] .= h(p, collect(Float64, -dt*(Ndelay-1):dt:(-dt)))
    end
    
    x0 = copy(x[Ndelay])
    
    if iip
        tmp1 = zeros(T, l)
        tmp2 = zeros(T, l)
        f(tmp1, x0, x[1], p, 0)
        f(tmp2, x[Ndelay+1], x[1], p, 0)
        @. x1[Ndelay+1] = x0 + dt^order*tmp1/(gamma(order)*order)
        @. x[Ndelay+1] = x0 + dt^order*(tmp2 + order*tmp1)/gamma(order+2)
    else
        @. x1[Ndelay+1] = x0 + dt^order*f(x0, x[1], p, 0)/(gamma(order)*order)
        @. x[Ndelay+1] = x0 + dt^order*(f(x[Ndelay+1], x[1], p, 0) + order*f(x0, x[1], p, 0))/gamma(order+2)
    end

    return ABMCache{iip, T}(prob, alg, mesh, u0, order, constant_lags, p, x, x0, x1, N, Ndelay, dt, kwargs)
    
end

function SciMLBase.solve!(cache::ABMCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, constant_algs, p, x, x0, x1, N, Ndelay, dt, kwargs = cache
    l = length(u0)
    if iip
        @fastmath @inbounds @simd for n=1:N
            tmp1 = zeros(T, l)
            prob.f(tmp1, x0, x[1], p, 0)
            M1 = @. (n^(order+1)-(n-order)*(n+1)^order)*tmp1
            N1 = @. ((n+1)^order-n^order)*tmp1

            @fastmath @inbounds @simd for j=1:n
                tmp2 = zeros(T, l)
                prob.f(tmp2, x[Ndelay+j], x[j], p, 0)
                M1 = @. M1+((n-j+2)^(order+1)+(n-j)^(order+1)-2*(n-j+1)^(order+1))*tmp2
                N1 = @. N1+((n-j+1)^order-(n-j)^order)*tmp2
            end
            @. x1[Ndelay+n+1] = x0+dt^order*N1/(gamma(order)*order)
            tmp3 = zeros(T, l)
            prob.f(tmp3, x[Ndelay+n+1], x[n+1], p, 0)
            @. x[Ndelay+n+1] = x0+dt^order*(tmp3 + M1)/gamma(order+2)
        end
    else
        @fastmath @inbounds @simd for n=1:N 
            @. M1=(n^(order+1)-(n-order)*(n+1)^order)*prob.f(x0, x[1], p, 0)
            @. N1=((n+1)^order-n^order)*prob.f(x0, x[1], p, 0)

            @fastmath @inbounds @simd for j=1:n
                @. M1 = M1+((n-j+2)^(order+1)+(n-j)^(order+1)-2*(n-j+1)^(order+1))*prob.f(x[Ndelay+j], x[j], p, 0)
                @. N1 = N1+((n-j+1)^order-(n-j)^order)*prob.f(x[Ndelay+j], x[j], p, 0)
            end
            @. x1[Ndelay+n+1] = x0+dt^order*N1/(gamma(order)*order)
            @. x[Ndelay+n+1] = x0+dt^order*(prob.f(x[Ndelay+n+1], x[n+1], p, 0)+M1)/gamma(order+2)
        end
    end
    
    u = x[Ndelay+1:end]

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end


#=
DelayABM method for system of fractional delay differential equations.

```tex
@inproceedings{Bhalekar2011APS,
  title={A PREDICTOR-CORRECTOR SCHEME FOR SOLVING NONLINEAR DELAY DIFFERENTIAL EQUATIONS OF FRACTIONAL ORDER},
  author={Sachin Bhalekar and Varsha Daftardar-Gejji},
  year={2011}
}
```
=#

function solve(FDDESys::FDDESystem, dt, ::DelayABM)
    @unpack f, ϕ, α, τ, T = FDDESys
    len = length(ϕ)
    N::Int = round(Int, T/dt)
    Ndelay = round(Int, τ/dt)
    t = collect(Float64, 0:dt:(T-τ))
    x1 = zeros(Ndelay+N+1, len)
    x = zeros(Ndelay+N+1, len)
    du = zeros(len)

    # Put the delay term in the array
    for i=1:len
        x[1:Ndelay, i] = ϕ[i]*ones(Ndelay)
    end
    
    x0 = copy(x[Ndelay, :])
    
    for i=1:len
        αi = α[i]
        f(du, x0, x[1, :], 0)
        x1[Ndelay+1, i] = x0[i] + dt^αi*du[i]/(gamma(αi)*αi)
    end

    for i=1:len
        αi = α[i]
        f(du, x[Ndelay+1, :], x[1, :], 0)
        du1 = copy(du)
        f(du, x0, x[1, :], 0)
        x[Ndelay+1, i] = x0[i] + dt^αi*(du1[i] + αi*du[i])/gamma(αi+2)
    end

    M1=zeros(len); N1=zeros(len)
    @fastmath @inbounds @simd for n=1:N
        for i=1:len
            αi = α[i]
            f(du, x0, x[1, :], 0)
            M1[i]=(n^(αi+1)-(n-αi)*(n+1)^αi)*du[i]
            N1[i]=((n+1)^αi-n^αi)*du[i]
        end

        @fastmath @inbounds @simd for j=1:n
            for i=1:len
                αi = α[i]
                f(du, x[Ndelay+j, :], x[j, :], 0)
                M1[i] = M1[i]+((n-j+2)^(αi+1)+(n-j)^(αi+1)-2*(n-j+1)^(αi+1))*du[i]
                N1[i] = N1[i]+((n-j+1)^αi-(n-j)^αi)*du[i]
            end
        end
        
        for i=1:len
            αi = α[i]
            x1[Ndelay+n+1, i] = x0[i]+dt^αi*N1[i]/(gamma(αi)*αi)
            f(du, x[Ndelay+n+1, :], x[n+1, :], 0)
            x[Ndelay+n+1, i] = x0[i]+dt^αi*(du[i]+M1[i])/gamma(αi+2)
        end
    end
    
    xresult = zeros(N-Ndelay+1, len)
    yresult = zeros(N-Ndelay+1, len)
    
    xresult[N-Ndelay+1, :] .= 0
    yresult[N-Ndelay+1, :] .= 0
    
    @fastmath @inbounds @simd for n=2*Ndelay+1:N+Ndelay+1  
       xresult[n-2*Ndelay, :] = x[n, :]
       yresult[n-2*Ndelay, :] = x[n-Ndelay, :]
    end

    #return xresult, yresult
    return FDDESystemSolution(t, yresult')
end