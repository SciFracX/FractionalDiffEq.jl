"""
# Usage

    solve(prob::FODESystem, h, GL())

Use Grunwald Letnikov difference method to solve system of system of FODE.

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
# Grunwald Letnikov discretization method dispatch for FODESystem
# struct GLWithMemory <: FractionalDiffEqAlgorithm end

function solve(prob::FODESystem, h, ::GL)
    # GL method is only for same order FODE
    @unpack f, α, u0, tspan, p = prob
    t0 = tspan[1]; T = tspan[2]
    α = α[1]
    hα = h^α[1]
    n::Int64 = floor(Int64, (T-t0)/h)+1
    l = length(u0)

    # Initialize solution
    result = zeros(Float64, length(u0), n)
    result[:, 1] = u0

    # generating generalized binomial Cα
    Cα = zeros(Float64, n)
    Cα[1] = 1

    @fastmath @inbounds @simd for j in range(2, n, step=1)
        Cα[j] = (1-(1+α)/(j-1))*Cα[j-1]
    end

    du = zeros(Float64, l)

    @fastmath @inbounds @simd for k in range(2, n, step=1)
        summation = zeros(Float64, length(u0))

        @fastmath @inbounds @simd for j in range(1, k-1, step=1)
            for i in eachindex(summation)
                summation[i] += Cα[j+1]*result[i, k-j]
            end
        end

        f(du, result[:, k-1], p, t0+(k-1)*h)
        result[:, k] = @. hα*du-summation
    end
    return result
end