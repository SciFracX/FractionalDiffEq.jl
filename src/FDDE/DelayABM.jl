"""
    solve(FDDE::FDDEProblem, h, DelayABM())

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
#FIXME: Also the problem definition f(t, ϕ, y) or f(t, y, ϕ)?
function solve(FDDE::FDDEProblem, h, ::DelayABM)
    @unpack f, ϕ, α, τ, tspan = FDDE
    N::Int = round(Int, tspan/h)
    Ndelay::Int = round(Int, τ/h)
    x1 = zeros(Float64, Ndelay+N+1)
    x = zeros(Float64, Ndelay+N+1)
    #x1[Ndelay+N+1] = 0
    
    #x[Ndelay+N+1] = 0

    # History function handling
    #FIXME: When the value of history function ϕ is different with the initial value?
    if typeof(ϕ) <: Number
        x[1:Ndelay] = ϕ*ones(Ndelay)
    elseif typeof(ϕ) <: Function
        x[Ndelay] = ϕ(0)
        x[1:Ndelay-1] .= ϕ(collect(Float64, -h*(Ndelay-1):h:(-h)))
    end
    
    x0 = copy(x[Ndelay])
    
    
    x1[Ndelay+1] = x0 + h^α*f(0, x[1], x0)/(gamma(α)*α)

    x[Ndelay+1] = x0 + h^α*(f(0, x[1], x[Ndelay+1]) + α*f(0, x[1], x0))/gamma(α+2)
    
    @fastmath @inbounds @simd for n=1:N 
        M1=(n^(α+1)-(n-α)*(n+1)^α)*f(0, x[1], x0)
        
        N1=((n+1)^α-n^α)*f(0, x[1], x0)

        @fastmath @inbounds @simd for j=1:n   
            M1 = M1+((n-j+2)^(α+1)+(n-j)^(α+1)-2*(n-j+1)^(α+1))*f(0, x[j], x[Ndelay+j])
            N1 = N1+((n-j+1)^α-(n-j)^α)*f(0, x[j], x[Ndelay+j])
        end
        x1[Ndelay+n+1] = x0+h^α*N1/(gamma(α)*α)
        x[Ndelay+n+1] = x0+h^α*(f(0, x[n+1], x[Ndelay+n+1])+M1)/gamma(α+2)
    end
    
    xresult = zeros(N-Ndelay+1)
    yresult = zeros(N-Ndelay+1)
    
    xresult[N-Ndelay+1]=0
    yresult[N-Ndelay+1]=0

    @fastmath @inbounds @simd for n=2*Ndelay+1:N+Ndelay+1  
       xresult[n-2*Ndelay] = x[n]
       yresult[n-2*Ndelay] = x[n-Ndelay]
    end

    return xresult, yresult
end