"""
# Usage

    solve(order, α, T, N, X, i, κ, v, fx0, fgz, f0t, flt, ADV_DIF())

### References

```tex
C.P. Li, R.F. Wu, H.F.  Ding. High-order approximation to Caputo derivatives and Caputo-type advection-diffusion equations (I). Communications in Applied and Industrial Mathematics, 2014, 6(2), e-536: 1-32.  DOI: 10.1685/journal.caim.536.
J.X. Cao, C.P. Li, Y.Q. Chen. High-order approximation to Caputo derivatives and Caputo-type advection-diffusion equations (II). Fractional Calculus and Applied Analysis, 2015, 18(3): 735-761.
H.F. Li, J.X. Cao, C.P. Li. High-order approximation to Caputo derivatives and Caputo-type advection-diffusion equations (III).submitted.
```
"""
struct ADV_DIF <: FPDEAlgorithm end

function solve(order::Integer, α, T, N, X, i, κ, v, fx0, fgz, f0t, flt, ::ADV_DIF)
    h = X/i
    gm = T/N
    miu = h^2*gm^(-α)/gamma(1-α)
    u = zeros(i-1, N+1)

    for i2=1:i-1
        u[i2, 1]=fx0(X*i2/i)
    end

    f = zeros(i-1, N)
    for i5=1:N
        for i6=1:i-1
            q1=[T*i5/N X*i6/i α]
            f[i6, i5]=fgz(q1)
        end
    end

    H = zeros(i-1, N)
    for i7 = 1:N
        q1 = [i7/N α]
        q2 = [T*(i7)/N α]
        H[1, i7] = (κ+v*h/2)*f0t(q1)
        H[i-1, i7] = (κ-v*h/2)*flt(q2)
    end

    for n=1:order-2
        g1 = ones(1, n)
        g2=g1*wj(n, n+1, α)
        g = zeros(n, 1)
        for i8=1:n
            g[i8, 1]=g2[1, i8]
        end
        u1=zeros(i-1, n)
        for i9=1:n
            u1[:, i9] = u[:, i9]
        end
        H1 = -miu*u1*g + h^2*f[:, n]+H[:, n]

        A = zeros(i-1, i-1)
        for i10=1:i-1
            A[i10, i10]=miu*g2[1, n+1]+2*κ
            if i10>1
                A[i10, i10-1]=-(κ+v*h/2)
            end
            if i10<i-1
                A[i10, i10+1]=-κ+v*h/2
            end
        end
        u[:, n+1]=zg(A, H1)
    end
    for n=order-1:N
        g1 = ones(1, n)
        g2 = g1*wj(n, order, α)
        g=zeros(n,1)
        for i8=1:n
            g[i8, 1]=g2[1, i8]
        end
        u1=zeros(i-1, n)
        for i9=1:n
            u1[:,i9] = u[:, i9]
        end
        H1=-miu*(u1*g) + h^2*f[:, n] + H[:, n]

        A=zeros(i-1, i-1)
        for i10=1:i-1
            A[i10, i10] = miu*g2[1, n+1] + 2*κ
            if i10 !== 1
                A[i10, i10-1] = -(κ+v*h/2)
            end
            if i10<i-1
                A[i10, i10+1]=-κ+v*h/2
            end
        end
        u[:, n+1] = zg(A, H1)
    end
    return u
end

function wj(n::Int, r::Int, a)
    A = zeros(n, n+1)
    for iw=1:r-2
        for jw=1:iw+1
        A[iw, jw] = w(iw+1-jw, iw+1, iw, n, a)
        end
    end
    for iw=r-1:n
        for jw=iw-r+2:iw+1
        A[iw, jw] = w(iw+1-jw, r, iw, n, a)
        end
    end
    return A
end

function zg(A, H)
    i = length(H)
    r = zeros(i)
    l = zeros(i)
    y = zeros(i)
    U = zeros(i)
    for t=1:i
        if t==1
            r[t] = A[1, 1]
        else
            l[t] = A[t, t-1]/r[t-1]
            r[t] = A[t, t] - l[t]*A[t-1, t]
        end
    end

    for p=1:i
        if p == 1
            y[p] = H[1, 1]
        else
            y[p] = H[p, 1] - l[p]*y[p-1]
        end
    end
    for k=0:i-1
       if k==0
        U[i] = y[i]/r[i]
       else
        U[i-k]=(y[i-k] - A[i-k, i-k+1]*U[i-k+1])/r[i-k]
       end
    end
    return U
end

function w(i::Int, r, j, n, a)
    ar=ones(r-1)
    br=ones(r-1)
    for lj=1:r-2
        jj=r-lj-1
        kj=r-2
        tj=i-1
        aj=collect(Int, 0:tj)
        bj=collect(Int, tj+2:kj+1)
        cj=[aj; bj]
        dj=collect(Int, -1:tj-1)
        ej=collect(Int, tj+1:kj)
        fj=[dj; ej]
        yj = binomial.(cj, jj)
        pj = binomial.(fj, jj)
        sj=1
        tj=1
        for m=1:jj
            sj = sj.*yj[:, m]
            tj = tj.*pj[:, m]
            ar[lj] = sum(sj)
            br[lj] = sum(tj)
        end
    end
    s=0
    for l=1:r-1
        k=1
        for m1=1:l
            k=k/(m1-a)
        end
        k = k*factorial(l)
        k = k*(ar[l]*((n-j)^(l-a))-br[l]*((n-j+1)^(l-a)))
        s = s+k
    end
    s = s*((-1)^(i+1))/(factorial(i)*factorial(r-1-i))
    return s
end