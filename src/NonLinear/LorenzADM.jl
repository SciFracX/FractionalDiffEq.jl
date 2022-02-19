struct FractionalLorenz
    a
    b
    c
    d
    α
end

struct LorenzADM <: FractionalDiffEqAlgorithm end

function solve(FL::FractionalLorenz, h, T, ::LorenzADM)
    a, b, c, d, q = FL.a, FL.b, FL.c, FL.d, FL.α
    N=round(Int, T/h)

    #k1, k2, k3 = zeros(n+1, 7), zeros(n+1, 7), zeros(n+1, 7)
    x1, x2, x3 = zeros(N+1), zeros(N+1), zeros(N+1)

    x1[1]=1
    x2[1]=2
    x3[1]=3

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
        k[2, 4]=c*k[1, 3]+d*k[2, 3]-k[1, 1]*k[3, 3]-k[1, 2]*k[3, 2]*gamma(2*q+1)/gamma(q+1)^2-k[1, 3]*k[3, 1]
        k[3, 4]=k[1, 1]*k[2, 3]+k[1, 2]*k[2, 2]*gamma(2*q+1)/gamma(q+1)^2+k[1, 3]*k[2, 1]-b*k[3, 3]

        k[1, 5]=a*(k[2, 4]-k[1, 4])
        k[2, 5]=c*k[1, 4]+d*k[2, 4]-k[1, 1]*k[3, 4]-k[1, 4]*k[3, 1]-(k[1, 3]*k[3, 2]+k[1, 2]*k[3, 3])*gamma(3*q+1)/(gamma(q+1)*gamma(2*q+1))
        k[3, 5]=k[1, 1]*k[2, 4]+k[1, 4]*k[2, 1]+b*k[3, 4]+(k[1, 3]*k[2, 2]+k[1, 2]*k[2, 3])*gamma(3*q+1)/(gamma(q+1)*gamma(2*q+1))

        k[1, 6]=a*(k[2, 5]-k[1, 5])
        k[2, 6]=c*k[1, 5]+d*k[2, 5]-k[1, 1]*k[3, 5]-(k[1, 4]*k[3, 2]+k[1, 2]*k[3, 4])*gamma(4*q+1)/(gamma(q+1)*gamma(3*q+1))-k[1, 3]*k[3, 3]*gamma(4*q+1)/gamma(2*q+1)^2-k[1, 5]*k[3, 1]
        k[3, 6]=k[1, 1]*k[2, 5]+(k[1, 4]*k[2, 2]+k[1, 2]*k[2, 4])*gamma(4*q+1)/(gamma(q+1)*gamma(3*q+1))+k[1, 3]*k[2, 3]*gamma(4*q+1)/gamma(2*q+1)^2+k[1, 5]*k[3, 1]

        k[1, 7]=a*(k[2, 6]-k[1, 6])
        k[2, 7]=c*k[1, 6]+d*k[2, 6]-k[1, 1]*k[3, 6]-(k[1, 2]*k[3, 5]+k[1, 5]*k[3, 2])*gamma(5*q+1)/(gamma(q+1)*gamma(4*q+1))-(k[1, 3]*k[3, 4]+k[1, 4]*k[3, 3])*gamma(5*q+1)/(gamma(2*q+1)*gamma(3*q+1))-k[1, 6]*k[3, 1]
        k[3, 7]=k[1, 1]*k[2, 6]+k[1, 6]*k[2, 1]-b*k[3, 6]+(k[1, 2]*k[2, 5]+k[1, 5]*k[2, 2])*gamma(5*q+1)/(gamma(q+1)*gamma(4*q+1))+(k[1, 3]*k[2, 4]+k[1, 4]*k[2, 3])*gamma(5*q+1)/(gamma(2*q+1)*gamma(3*q+1))


        temp1, temp2, temp3 =0, 0, 0
        for j=1:7
            temp1 += k[1, j]*h^((j-1)*q)/gamma((j-1)*q+1)
            temp2 += k[2, j]*h^((j-1)*q)/gamma((j-1)*q+1)
            temp3 += k[3, j]*h^((j-1)*q)/gamma((j-1)*q+1)
        end
        x1[n+1]=temp1
        x2[n+1]=temp2
        x3[n+1]=temp3
    end
    return x1, x2, x3
end