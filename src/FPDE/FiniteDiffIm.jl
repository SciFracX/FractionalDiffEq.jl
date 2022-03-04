"""
    solve(α, dx, dt, xStart, xEnd, n, κ, FiniteDiffIm())

Use implicit finite difference to discrete equations.

### References

```tex
@article{Murio2008ImplicitFD,
  title={Implicit finite difference approximation for time fractional diffusion equations},
  author={Diego A. Murio},
  journal={Comput. Math. Appl.},
  year={2008},
  volume={56},
  pages={1138-1145}
}
```

Matlab version: https://github.com/awstown/Fractional-Derivative
"""
struct FiniteDiffIm <: FractionalDiffEqAlgorithm end


function solve(α, dx, dt, xStart, xEnd, n, κ, u0t, uendt, u0, ::FiniteDiffIm)
    x = collect(0:dx:xEnd)
    t = collect(0:dt:n)

    U = zeros(Int64(n/dt + 1), Int64((xEnd - xStart)/dx + 1))

    # Boundry conditions
    U[:, 1] .= u0t
    U[:, end] .= uendt

    U[1, :] .= u0.(x)

    mu = (dt^α)/(dx^2)
    r = mu * gamma(2-α)
    A1 = diagm((1+2*r)*ones(length(x)-2))
    A2 = diagm(1 => -r*ones(length(x)-3))
    A3 = diagm(-1 => -r*ones(length(x)-3))
    A = A1 + A2 + A3


    U[2, 2:end-1] = A \ U[1, 2:end-1]

    j = collect(0:length(t))
    b_j = (j.+1).^(1-α) .- j.^(1-α)
    c_j = b_j[1:end-1] - b_j[2:end]

    V = copy(U[2:end, 2:end-1])
    for k = 1:length(V[:, 1])-1
        V[k+1, :] = A \ ImNextStep(V, U[1, 2:end-1], k, b_j, c_j)'
    end
    U[2:end, 2:end-1] = V
    return U
end

ImNextStep( V, U_0, k, b_j, c_j ) = c_j[1:k]' * reverse(V[1:k, :], dims=1) .+ (b_j[k+1].*U_0)'