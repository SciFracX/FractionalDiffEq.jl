using SpecialFunctions
"""
  solve(FDDE::FDDESystem, h, DelayABMYuan())

Solve system of fractional delay differential equations.

### References

```tex
@article{Yuan2013ChaosDA,
  title={Chaos detection and parameter identification in fractional-order chaotic systems with delay},
  author={Liguo Yuan and Qigui Yang and Caibin Zeng},
  journal={Nonlinear Dynamics},
  year={2013},
  volume={73},
  pages={439-448}
}
```
"""
struct DelayABMYuan <: FractionalDiffEqAlgorithm end

function solve(FDDESys::FDDESystem, h, ::DelayABMYuan)
  @unpack f, ϕ, α, τ, T = FDDESys
  N = round(Int, T/h)
  Ndelay = round(Int, τ/h)
  m = length(ϕ)

  x1=zeros(Ndelay+N+1, m)
  x1[Ndelay+N+1, :] .= 0 #Initialization

  x=zeros(Ndelay+N+1, m)
  x[Ndelay+N+1, :] .= 0

  for i=1:m
      x[1:Ndelay, i] = ϕ[i]*ones(Ndelay)
  end

  x0=copy(x[Ndelay, :])

  for i=1:m
    αi = α[i]
    x1[Ndelay+1, i]=x0[i]+h^αi*(f(x0..., x[1, :]..., i))/(gamma(αi)*αi)
  end

  for i=1:m
    αi = α[i]
    x[Ndelay+1, i]=x0[i]+h^αi*(f(x1[Ndelay+1, :]..., x[2, :]..., i)+αi*f(x0..., x[1, :]..., i))/gamma(αi+2)
  end

  M1=zeros(m)
  N1=zeros(m)
  for n=1:N
    
    for i=1:m
      αi = α[i]
      M1[i]=(n^(αi+1)-(n-αi)*(n+1)^αi)*(f(x0..., x[1, :]..., i))
      N1[i]=((n+1)^αi-n^αi)*(f(x0..., x[1, :]..., i))
    end

    for j=1:n
      for i=1:m
        αi = α[i]
        M1[i]=M1[i]+((n-j+2)^(αi+1)+(n-j)^(αi+1)-2*(n-j+1)^(αi+1))*(f(x[Ndelay+j, :]..., x[j, :]..., i))
        N1[i]=N1[i]+((n-j+1)^αi-(n-j)^αi)*(f(x[Ndelay+j, :]..., x[j, :]..., i))
      end
    end

    @. x1[Ndelay+n+1, :]=x0+h^α*N1/(gamma(α)*α)
  
    for i=1:m
      αi = α[i]
      x[Ndelay+n+1, i]=x0[i]+h^αi*(f(x[Ndelay+n+1, :]..., x[n+1, :]..., i)+M1[i])/gamma(αi+2)
    end
  end
  return x
end

