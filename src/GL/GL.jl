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

function solve(FODE::SingleTermFODEProblem, h, ::GL)
    @unpack f, α, u0, T = FODE
    N = floor(Int, T/h)+1
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
        right = 0
        for j=1:i-1
            right += c[j]*y[i-j]
        end
        y[i] = f((i-1)*h, y[i-1])*h^α - right
    end
    return y
end