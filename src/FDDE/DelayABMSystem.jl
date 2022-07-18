
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

function solve(FDDESys::FDDESystem, h, ::DelayABM)
    @unpack f, ϕ, α, τ, T = FDDESys
    len = length(ϕ)
    N::Int = round(Int, T/h)
    Ndelay = round(Int, τ/h)
    t = collect(Float64, 0:h:(T-τ))
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
        x1[Ndelay+1, i] = x0[i] + h^αi*du[i]/(gamma(αi)*αi)
    end

    for i=1:len
        αi = α[i]
        f(du, x[Ndelay+1, :], x[1, :], 0)
        du1 = copy(du)
        f(du, x0, x[1, :], 0)
        x[Ndelay+1, i] = x0[i] + h^αi*(du1[i] + αi*du[i])/gamma(αi+2)
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
            x1[Ndelay+n+1, i] = x0[i]+h^αi*N1[i]/(gamma(αi)*αi)
            f(du, x[Ndelay+n+1, :], x[n+1, :], 0)
            x[Ndelay+n+1, i] = x0[i]+h^αi*(du[i]+M1[i])/gamma(αi+2)
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