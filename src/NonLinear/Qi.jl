#import FractionalDiffEq: FractionalDiffEqAlgorithm, solve
#=
Example: https://ieeexplore.ieee.org/abstract/document/5566025/
=##=
"""
# Qi System

```math

```

```tex
@INPROCEEDINGS{5566025,
  author={Wu, Xiangjun and Yang, Yang}, 
  booktitle={2010 International Conference on Intelligent Computing and Cognitive Informatics},
  title={Chaos in the Fractional-Order Qi System and Its Synchronization Using Active Control},  
  year={2010},  volume={},  number={},  pages={109-112},  keywords={},  
  doi={10.1109/ICICCI.2010.63},  
  ISSN={},  
  month={June}
  }
```
"""
struct QiSystem
    parameters
    orders
    T
    Y0
end

struct QiAlg <: FractionalDiffEqAlgorithm end
=#
function solve(parameters, orders, T, Y0, h)
    n=Int64(floor(T/h))
    q1=orders[1]
    q2=orders[2]
    q3=orders[3]

    a=parameters[1]
    b=parameters[2]
    c=parameters[3]
    d=parameters[4]
    r=parameters[5]
    cp1=1
    cp2=1
    cp3=1
    c1=zeros(n)
    c2=zeros(n)
    c3=zeros(n)


    for j=1:n
        c1[j]=(1-(1+q1)/j)*cp1
        c2[j]=(1-(1+q2)/j)*cp2
        c3[j]=(1-(1+q3)/j)*cp3

        cp1=c1[j]
        cp2=c2[j]
        cp3=c3[j]
    end

    x=zeros(n)
    y=zeros(n)
    z=zeros(n)

    x[1]=Y0[1]
    y[1]=Y0[2]
    z[1]=Y0[3]



    for i=2:n
        x[i]=(-a*x[i-1]+a*y[i-1]+r*y[i-1]*z[i-1])*h^q1 - memo(x, c1, i)
        y[i]=(c*x[i-1]+d*y[i-1]-x[i-1]*z[i-1])*h^q2 - memo(y, c2, i)
        z[i]=(-b*z[i-1]+x[i-1]*y[i-1])*h^q3 - memo(z, c3, i)
    end

    return x, y, z
end
function memo(r, c, k)
    temp=1
    for j=1:k-1
        temp = temp+c[j]*r[k-j]
    end
    return temp
end

#=

x, y, z = testsolve([35, 8/3, 80, -1, 1], [0.96, 0.96, 0.96], 50, [0.1, 0.2, 0.3], 0.0001)
using Plots

plot(x, y)
=#