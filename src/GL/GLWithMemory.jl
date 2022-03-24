"""
# Usage

    solve(prob::FODESystem, h, T, GLWithMemory())

### Reference

```tex
@INPROCEEDINGS{8742063,  
author={Clemente-López, D. and Muñoz-Pacheco, J. M. and Félix-Beltrán, O. G. and Volos, C.},  
booktitle={2019 8th International Conference on Modern Circuits and Systems Technologies (MOCAST)},   
title={Efficient Computation of the Grünwald-Letnikov Method for ARM-Based Implementations of Fractional-Order Chaotic Systems},   
year={2019},   
doi={10.1109/MOCAST.2019.8742063}}
```

Python version by https://github.com/DClementeL/Grunwald_Letnikov
"""
struct GLWithMemory <: FractionalDiffEqAlgorithm end

function solve(prob::FODESystem, h, tf, ::GLWithMemory)
    @unpack f, α, x0 = prob
    hα=h^α[1]
    n::Int64 = floor(Int64, tf/h)+1
    l = length(x0)

    # Initialize solution
    result = zeros(n, length(x0))
    result[1, :] = x0

    # generating coefficients Cα
    Cα = zeros(n)
    Cα[1] = 1
    @fastmath @inbounds @simd for j in range(2, n, step=1)
        Cα[j] = (1-(1+α[1])/(j-1))*Cα[j-1]
    end

    du = zeros(l)

    @fastmath @inbounds @simd for k in range(2, n, step=1)
        summation = zeros(length(x0))

        @fastmath @inbounds @simd for j in range(1, k-1, step=1)
            for i in eachindex(summation)
                summation[i] += Cα[j+1]*result[k-j, i]
            end
        end

        f(du, result[k-1, :], nothing, (k-1)*h)
        for i in range(1, l, step=1)
            result[k, i] = hα*du[i] - summation[i]
        end
    end
    return result
end