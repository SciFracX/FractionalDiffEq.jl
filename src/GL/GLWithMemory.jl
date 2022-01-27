import FractionalDiffEq: FractionalDiffEqAlgorithm, solve

"""
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

function solve(f, α, x0, h, tf, ::GLWithMemory)

    hα=h^α[1]

    n = Int64(tf/h+1)
    x, y, z = zeros(n), zeros(n), zeros(n)
    x[1], y[1], z[1]= x0[1], x0[2], x0[3]

    x1 = copy(x0)

    Cα = zeros(n)
    Cα[1] = 1

    for j in range(2, n, step=1)
        Cα[j] = (1-(1+α[1])/(j-1))*Cα[j-1]
    end

    for k in range(2, n, step=1)
        sum1, sum2, sum3 = 0, 0, 0

        for j in range(1, k-1, step=1)
            sum1 += Cα[j+1]*x[k-j]
            sum2 += Cα[j+1]*y[k-j]
            sum3 += Cα[j+1]*z[k-j]
        end

            x[k]=hα*f((k-1)*h, x[k-1], y[k-1], z[k-1], 1) - sum1
            y[k]=hα*f((k-1)*h, x[k-1], y[k-1], z[k-1], 2) - sum2
            z[k]=hα*f((k-1)*h, x[k-1], y[k-1], z[k-1], 3) - sum3
    end
    return x, y, z
end