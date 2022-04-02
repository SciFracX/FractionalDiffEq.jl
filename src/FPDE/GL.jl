"""
Using Grunwald Letnikov discretization operator to discrete fractional operator.
"""
struct GLDiff <: FractionalDiffEqAlgorithm end

function solve(α, d, rightfun, M, N, initial_condition, left_boundry, right_boundry, ::GLDiff)
    ue=zeros(M+1, N+1)
    u=copy(ue)
    D=zeros(M-1, 1)
    a=copy(D)
    #

    #exact(x, t) = exp(-t).*x.^3
    #=
    for k=1:N+1
        ue[1:end, k]=exact(x[1:end], t[k])
    end
    =#
    u[1:end, 1]=initial_condition.(x)
    u[1, 1:end]=left_boundry.(t)
    u[end, 1:end]=right_boundry.(t)
    A=zeros(M-1, M-1)
    for i=1:M-1
        D[i, 1]=d(x[i+1])
    end
    a=tau*D/(2*h^α)
    gg=g(M, α)
    for i=1:M-1
        for k=1:N-1
            if k <= i-1
                A[i, k]=a[i, 1]*gg[i-k+2, 1]
            elseif k==i
                A[i, k]=a[i, 1]*gg[2, 1]
            elseif k==i+1
                A[i, k]=a[i, 1]*gg[1, 1]
            else
                A[i, k]=0
            end
        end
    end

    for k=1:N
        b = (zeros(M-1, M-1)+I+A)*u[2:end-1, k].+tau.*rightfun.(x[2:end-1], t[k]+tau/2) .+a.*(gg[3:end].*(u[1, k+1]+u[1, k]).*ones(M-1, 1))
        b[end, 1]=b[end, 1]+a[end, 1]*gg[1]*(u[end, k+1]+u[end, k])
        u[2:end-1, k+1]=(zeros(M-1, M-1)+I-A)\b
    end
    return u
end

function g(M, α)
    gg_alpha = zeros(M+1, 1)
    gg_alpha[1, 1]=1
    for i=1:M
        gg_alpha[i+1, 1]=gamma(i-α)/(gamma(-α)*gamma(i+1))
    end
    return gg_alpha
end
