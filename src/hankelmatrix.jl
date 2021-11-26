using SpecialMatrices

function solve(a, na, b, nb, u, t)
    h=t[2]-t[1]
    u=u[:]
    A, B = 0, 0

    g=genfun(1)
    nt=length(t)
    n=length(a)
    m=length(b)
    for i=1:n
        A=A+getvec(na[i], nt, g)*a[i]/(h^na[i])
    end

    for i=1:m
        B=B+getvec(nb[i], nt, g)*b[i]/(h^nb[i])
    end

    A=rotl90(newhankel(A[end:-1:1]))
    B=rotl90(newhankel(B[end:-1:1]))

    y=B*inv(A)*u
    return y
end


"""
P-th precision polynomial generate function
```math
g_p(z)=\\sum_{k=1}^p \\frac{1}{k}(1-z)^k
```
"""
function genfun(p)
    a=collect(1:p+1)
    A=Vandermonde(a)'
    return (1 .-a')*inv(A')
end

function getvec(α, n, g)
    p=length(g)-1
    b=1+α
    g0=g[1]
    w=zeros(1, n)
    w[1]=g[1]^α

    for m=2:p
        M=m-1
        dA=b/M
        w[m]=(-(g[2:m].*collect((1-dA):-dA:(1-b)))*w[M:-1:1]'/g0)[1]
    end

    for k=p+1:n
        M=k-1
        dA=b/M
        w[k]=(-(g[2:(p+1)].*collect((1-dA):-dA:(1-p*dA)))*w[M:-1:(k-p)]'/g0)[1]
    end

    return w
end

function newhankel(v)
    n=length(v)
    v=v[:]

    hankelm=zeros(n, n)
    for i=1:length(v)
        hankelm[i, 1:end-i+1]=v[i:end]
    end

    return hankelm

end


t=collect(0:0.002:10);

u=sin.(t.^2);

result=solve([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], [30 90], [1 0.3], u, t);

plot(t, result)


# It need to known that in Julia, for is somewhat slowerer than vectorization