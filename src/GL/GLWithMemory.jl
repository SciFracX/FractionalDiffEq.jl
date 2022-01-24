import FractionalDiffEq.FractionalDiffEqAlgorithm


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
#=
function solve()
    α=0.99
    h=0.005
    hα=h^α
    tf=30
    n=Int64(tf/h+1)
    x, y, z = zeros(n), zeros(n), zeros(n)
    x[1]=1
    y[1]=0
    z[1]=1

    a, b, c = 10, 28, 8/3

    Cα=zeros(n)
    Cα[1]=1

    for j in range(2, n, step=1)
        Cα[j] = (1-(1+α)/(j-1))*Cα[j-1]
    end

    for k in range(2, n, step=1)
        sum1, sum2, sum3 = 0, 0, 0

        for j in range(1, k-1, step=1)
            sum1 += Cq[j+1]*x[k-j]
            sum2 += Cq[j+1]*y[k-j]
            sum3 += Cq[j+1]*z[k-j]
        end

        x[k]=hα*(a*(y[k-1]-x[k-1]))-sum1
        y[k]=hα*(x[k-1]*(b-z[k-1]) - y[k-1])-sum2
        z[k]=hα*(x[k-1]*y[k-1]-c*z[k-1])-sum3
    end
    return x, y, z
end
=#