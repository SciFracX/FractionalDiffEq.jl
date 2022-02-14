"""
    solve()

Use explicit Caputo discretization method

### References

```tex
@article{Murillo2011AnED,
  title={An Explicit Difference Method for Solving Fractional Diffusion and Diffusion-Wave Equations in the Caputo Form},
  author={Joaqu{\'i}n Quintana Murillo and Santos B. Yuste},
  journal={Journal of Computational and Nonlinear Dynamics},
  year={2011},
  volume={6},
  pages={021014}
}
```

Matlab version: https://github.com/awstown/Fractional-Derivative
"""
struct CaputoDiscretizationEX <: FractionalDiffEqAlgorithm end

function solve(fdorder, dx, dt, xStart, xEnd, n, K, ::CaputoDiscretizationEX)
    x = collect(0:dx:xEnd)
    t = collect(0:dt:n)
    S = K * ((dt^fdorder)/(dx^2))
    S_bar = gamma(3 - fdorder) * S
        

    U = zeros(Int64(n/dt + 1), round(Int, (xEnd - xStart)/dx + 1))

    # Boundry conditions handling
    U[:, 1] .= 0
    U[:, end] .= 0
    U[1, :] = sin.(x)

    k = collect(1:length(t)-2)
    bOfK = bbcoeff(k, fdorder)
    test = Float64[]

    for m = 1:(length(t)- 1)
        # It is from 2 to end - 1 because of the B.C.
        U[m+1, 2:end-1], diff = nextStep(U, m, S_bar, test, bOfK);
    end
    return U
end



function bbcoeff(k, fdorder)
    c = zeros(length(k))
    for i=1:length(k)
        if maximum(k) < 2
            @. c = (k+1)^(1-fdorder) - k^(1-fdorder);
        else
            @. c = (k+1)^(2-fdorder) - k^(2-fdorder);
        end
    end
    return c
end

# Compute values for each time step
function nextStep(U, m, S_bar, diff, bOfK )
    row_left   = U[m, 1:end-2]
    row_right  = U[m, 3:end]
    row_center = U[m, 2:end-1]
    
    if m < 2
        Uspatial = row_center .+ S_bar.*(row_left - 2*row_center + row_right)
    else
        row_below  = U[m-1, 2:end-1]
        Uspatial = 2 .*row_center .- row_below .+ S_bar.*(row_left .- 2*row_center .+ row_right);
    end

    if m > 2
        if length(diff) == 0
            diff = U[m, 2:end-1] - 2*U[m-1, 2:end-1] + U[m-2, 2:end-1]
        else
            diff = [U[m, 2:end-1] .- 2 .*U[m-1, 2:end-1] .+ U[m-2, 2:end-1]; diff]
        end
        Utemporal = bOfK[1:length(diff[:, 1])]' * diff
        Unext = Uspatial .- Utemporal;
    elseif m > 1
        if length(diff) == 0
            diff = U[m, 2:end-1] .- U[m-1, 2:end-1]
        else
            diff = [U[m, 2:end-1] .- U[m-1,2:end-1]; diff]
        end
        Utemporal = bOfK[1:length(diff[:, 1])]' * diff
        Unext = Uspatial .- Utemporal
    else
        Unext = Uspatial
    end

    return Unext, diff
end