using SpecialMatrices

function nlsolve(f, α, x0, h, tn)    
    n = length(x0)
    m = Int64(round(tn/h)+1)
    g=genfun(1)
    x0=x0[:]
    ha=h.^α
    z=zeros(n, m)
    x1=x0

    W = zeros(n, m)

    for i=1:n
        W[i, :]=getvec(α[i], m, g)
    end

    for k=2:m
        tk=(k-1)*h
        L=k-1
        for i=1:n
            x1[i]=f(tk, x1, i)*ha[i] - W[i, 2:L+1]'*z[i, k-1:-1:k-L] + x0[i]
        end
        z[:, k] = x1-x0
    end

    x=(z + repeat(x0, 1, m))'
    return x

end

function genfun(p)
    a = collect(1:p+1)
    A = Vandermonde(a)'
    return (1 .-a')*inv(A')
end

function getvec(α, n, g)
    p = length(g)-1
    b = 1 + α
    g0 = g[1]
    w = Float64[]
    push!(w, g[1]^α)

    for m = 2:p
        M = m-1
        dA = b/M
        temp = (-(g[2:m] .*collect((1-dA):-dA:(1-b))))' *w[M:-1:1]/g0
        push!(w, temp)
    end

    for k = p+1:n
        M = k-1
        dA = b/M
        temp = (-(g[2:(p+1)] .*collect((1-dA):-dA:(1-p*dA))))' *w[M:-1:(k-p)]/g0
        push!(w, temp)
    end

    return w
end

function chua(t, x, k)
    a=10.725
    b=10.593
    c=0.268
    m0=-1.1726
    m1=-0.7872

    if k==1
        f=m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))
        y=a*(x[2]-x[1]-f)
        return y
    elseif k==2
        y=x[1]-x[2]+x[3]
        return y
    elseif k==3
        y=-b*x[2]-c*x[3]
        return y
    end
end

alpha=[0.93, 0.99, 0.92]
x0=[0.2; -0.1; 0.1]
h=0.1
tn=2
result=nlsolve(chua, alpha, x0, h, tn)

#print(result)


m=result[1, :]
n=result[2, :]

using Plots
gr()
plot(n, m)
