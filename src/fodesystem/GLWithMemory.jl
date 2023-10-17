"""
# Usage

    solve(prob::FODEProblem, h, GL())

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
# Grunwald Letnikov discretization method dispatch for FODEProblem
# struct GLWithMemory <: FractionalDiffEqAlgorithm end
struct GL <: FODESystemAlgorithm end

function solve(prob::FODEProblem, h, ::GL)
    # GL method is only for same order FODE
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; T = tspan[2]
    t = collect(Float64, t0:h:T)
    order = order[1]
    horder = h^order[1]
    n::Int64 = floor(Int64, (T-t0)/h)+1
    l = length(u0)

    # Initialize solution
    result = zeros(Float64, length(u0), n)
    result[:, 1] = u0

    # generating generalized binomial Corder
    Corder = zeros(Float64, n)
    Corder[1] = 1

    @fastmath @inbounds @simd for j in range(2, n, step=1)
        Corder[j] = (1-(1+order)/(j-1))*Corder[j-1]
    end

    du = zeros(Float64, l)

    @fastmath @inbounds @simd for k in range(2, n, step=1)
        summation = zeros(Float64, length(u0))

        @fastmath @inbounds @simd for j in range(1, k-1, step=1)
            for i in eachindex(summation)
                summation[i] += Corder[j+1]*result[i, k-j]
            end
        end

        f(du, result[:, k-1], p, t0+(k-1)*h)
        result[:, k] = @. horder*du-summation
    end
    return FODESystemSolution(t, result)
end