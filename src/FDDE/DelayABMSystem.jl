
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
    x1 = zeros(Ndelay+N+1, len)
    x = zeros(Ndelay+N+1, len)
    #x1[Ndelay+N+1] = 0
    
    #x[Ndelay+N+1] = 0

    # History function handling
    #FIXME: When the value of history function ϕ is different with the initial value?
    for i=1:len
        x[1:Ndelay, i] = ϕ[i]*ones(Ndelay)
    end
    
    x0 = copy(x[Ndelay, :])
    
    for i=1:len
        x1[Ndelay+1, i] =  x0[i] + h^α*f(0, x[1, :], x0, i)/(gamma(α)*α)
    end

    for i=1:len
        x[Ndelay+1, i] = x0[i] + h^α*(f(0, x[1, :], x[Ndelay+1, :], i) + α*f(0, x[1, :], x0, i))/gamma(α+2)
    end
    M1=zeros(len); N1=zeros(len)
    @fastmath @inbounds @simd for n=1:N
        for i=1:len
            M1[i]=(n^(α+1)-(n-α)*(n+1)^α)*f(0, x[1, :], x0, i)
        end
        
        for i=1:len
            N1[i]=((n+1)^α-n^α)*f(0, x[1, :], x0, i)
        end

        @fastmath @inbounds @simd for j=1:n
            for i=1:len
                M1[i] = M1[i]+((n-j+2)^(α+1)+(n-j)^(α+1)-2*(n-j+1)^(α+1))*f(0, x[j, :], x[Ndelay+j, :], i)
                N1[i] = N1[i]+((n-j+1)^α-(n-j)^α)*f(0, x[j, :], x[Ndelay+j, :], i)
            end
        end
        
        for i=1:len
            x1[Ndelay+n+1, i] = x0[i]+h^α*N1[i]/(gamma(α)*α)
            x[Ndelay+n+1, i] = x0[i]+h^α*(f(0, x[n+1, :], x[Ndelay+n+1, :], i)+M1[i])/gamma(α+2)
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

    return xresult, yresult
end