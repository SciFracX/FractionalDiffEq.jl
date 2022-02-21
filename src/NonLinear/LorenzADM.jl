struct FractionalLorenz
    a
    b
    c
    d
    α
    x0
end

"""

### References

```tex
@article{He2017DynamicsOT,
  title={Dynamics of the Fractional-order Lorenz System Based on Adomian Decomposition Method and Its DSP Implementation},
  author={Shaobo He and Kehui Sun and Huihai Wang},
  journal={IEEE/CAA Journal of Automatica Sinica},
  year={2017},
  pages={0-0}
}
```
"""
struct LorenzADM <: FractionalDiffEqAlgorithm end

function solve(FL::FractionalLorenz, h, T, ::LorenzADM)
    @unpack a, b, c, d, α, x0 = FL
    N=round(Int, T/h)

    #k1, k2, k3 = zeros(n+1, 7), zeros(n+1, 7), zeros(n+1, 7)
    x1, x2, x3 = zeros(N+1), zeros(N+1), zeros(N+1)

    x1[1], x2[1], x3[1] = x0

    k=zeros(3, 7)
    for n=1:N
        k[:, 1] = [x1[n], x2[n], x3[n]]
        k[1, 2] = a*(k[2, 1]-k[1, 1])
        k[2, 2] = c*k[1, 1]+d*k[2, 1]-k[1, 1]*k[3, 1]
        k[3, 2] = -b*k[3, 1]+k[1, 1]*k[2, 1]

        k[1, 3]=a*(k[2, 2]-k[1, 2])
        k[2, 3]=c*k[1, 2]+d*k[2, 2]-k[1, 1]*k[3, 2]-k[1, 2]*k[3, 1]
        k[3, 3]=k[1, 2]*k[2, 1]+k[1, 1]*k[2, 2]-b*k[3, 2]

        k[1, 4]=a*(k[2, 1]-k[1, 1])
        k[2, 4]=c*k[1, 3]+d*k[2, 3]-k[1, 1]*k[3, 3]-k[1, 2]*k[3, 2]*gamma(2*α+1)/gamma(α+1)^2-k[1, 3]*k[3, 1]
        k[3, 4]=k[1, 1]*k[2, 3]+k[1, 2]*k[2, 2]*gamma(2*α+1)/gamma(α+1)^2+k[1, 3]*k[2, 1]-b*k[3, 3]

        k[1, 5]=a*(k[2, 4]-k[1, 4])
        k[2, 5]=c*k[1, 4]+d*k[2, 4]-k[1, 1]*k[3, 4]-k[1, 4]*k[3, 1]-(k[1, 3]*k[3, 2]+k[1, 2]*k[3, 3])*gamma(3*α+1)/(gamma(α+1)*gamma(2*α+1))
        k[3, 5]=k[1, 1]*k[2, 4]+k[1, 4]*k[2, 1]+b*k[3, 4]+(k[1, 3]*k[2, 2]+k[1, 2]*k[2, 3])*gamma(3*α+1)/(gamma(α+1)*gamma(2*α+1))

        k[1, 6]=a*(k[2, 5]-k[1, 5])
        k[2, 6]=c*k[1, 5]+d*k[2, 5]-k[1, 1]*k[3, 5]-(k[1, 4]*k[3, 2]+k[1, 2]*k[3, 4])*gamma(4*α+1)/(gamma(α+1)*gamma(3*α+1))-k[1, 3]*k[3, 3]*gamma(4*α+1)/gamma(2*α+1)^2-k[1, 5]*k[3, 1]
        k[3, 6]=k[1, 1]*k[2, 5]+(k[1, 4]*k[2, 2]+k[1, 2]*k[2, 4])*gamma(4*α+1)/(gamma(α+1)*gamma(3*α+1))+k[1, 3]*k[2, 3]*gamma(4*α+1)/gamma(2*α+1)^2+k[1, 5]*k[3, 1]

        k[1, 7]=a*(k[2, 6]-k[1, 6])
        k[2, 7]=c*k[1, 6]+d*k[2, 6]-k[1, 1]*k[3, 6]-(k[1, 2]*k[3, 5]+k[1, 5]*k[3, 2])*gamma(5*α+1)/(gamma(α+1)*gamma(4*α+1))-(k[1, 3]*k[3, 4]+k[1, 4]*k[3, 3])*gamma(5*α+1)/(gamma(2*α+1)*gamma(3*α+1))-k[1, 6]*k[3, 1]
        k[3, 7]=k[1, 1]*k[2, 6]+k[1, 6]*k[2, 1]-b*k[3, 6]+(k[1, 2]*k[2, 5]+k[1, 5]*k[2, 2])*gamma(5*α+1)/(gamma(α+1)*gamma(4*α+1))+(k[1, 3]*k[2, 4]+k[1, 4]*k[2, 3])*gamma(5*α+1)/(gamma(2*α+1)*gamma(3*α+1))


        temp1, temp2, temp3 =0, 0, 0
        for j=1:7
            temp1 += k[1, j]*h^((j-1)*α)/gamma((j-1)*α+1)
            temp2 += k[2, j]*h^((j-1)*α)/gamma((j-1)*α+1)
            temp3 += k[3, j]*h^((j-1)*α)/gamma((j-1)*α+1)
        end
        x1[n+1]=temp1
        x2[n+1]=temp2
        x3[n+1]=temp3
    end
    return x1, x2, x3
end