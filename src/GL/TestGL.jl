import FractionalDiffEq: FractionalDiffEqAlgorithm, SingleTermFODEProblem, solve

"""
Grunwald Letnikov method for fractional ordinary differential equations

```tex
@INPROCEEDINGS{8742063,  
author={Clemente-López, D. and Muñoz-Pacheco, J. M. and Félix-Beltrán, O. G. and Volos, C.},  
booktitle={2019 8th International Conference on Modern Circuits and Systems Technologies (MOCAST)},   
title={Efficient Computation of the Grünwald-Letnikov Method for ARM-Based Implementations of Fractional-Order Chaotic Systems},   
year={2019},   
doi={10.1109/MOCAST.2019.8742063}}
```
"""
struct GL <: FractionalDiffEqAlgorithm end

# Some points are a little big
function solve(prob::SingleTermFODEProblem, u0, T, ::GL)
    f, α, h = prob.f, prob.α, prob.h
    N = Int64(floor(T/h))

    c = zeros(N)

    cp = 1
    for j = 1:N
        c[j] = (1-(1+α)/j)*cp
        cp = c[j]
    end

    # Initialization
    y = zeros(N)
    y[1] = u0

    for i = 2:N
        y[i] = f(y[i-1])*h^α - Cq(y, c, i)
    end

    return y
end

function Cq(r, c, k)
    temp = 1
    for j = 1:k-1
        temp = temp + c[j]*r[k-j]
    end
    return temp
end