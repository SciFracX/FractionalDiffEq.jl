struct ADV_DIF <: FractionalDiffEqAlgorithm end

function solve(q, fx0, fgz, f0t, flt, ::ADV_DIF)
    r=Int64(q[1,1]);a=q[1,2];T=q[1,3];N=Int64(q[1,4]);X=q[1,5];i=Int64(q[1,6]);k=q[1,7];v=q[1,8]
    h=X/i;gm=T/N;miu=h^2*gm^(-a)/gamma(1-a)
    u = zeros(i-1, N+1)

    for i2=1:i-1
        u[i2, 1]=fx0(X*(i2)/i)
    end

    f = zeros(i-1, N)
    for i5=1:N
        for i6=1:i-1
            q1=[T*(i5)/N X*(i6)/i a]
            f[i6, i5]=fgz(q1)
        end
    end

    H=zeros(i-1,N)
    for i7=1:N
        q1 = [i7/N a]
        q2 = [T*(i7)/N a]
        H[1, i7] = (k+v*h/2)*f0t(q1)
        H[i-1, i7] = (k-v*h/2)*flt(q2)
    end

    for n=1:r-2
        gp=[n n+1 a]
        g1 = ones(1, n)
        g2=g1*wj(gp)
        g = zeros(n, 1)
        for i8=1:n
            g[i8, 1]=g2[1,i8]
        end
        u1=zeros(i-1, n)
        for i9=1:n
            u1[:, i9]=u[:, i9]
        end
        H1=-miu*(u1*g) + h^2*f[:, n]+H[:, n]

        A=zeros(i-1,i-1)
        for i10=1:i-1
            A[i10, i10]=miu*g2[1, n+1]+2*k
            if i10>1
                A[i10, i10-1]=-(k+v*h/2)
            end
           if i10<i-1
                A[i10, i10+1]=-k+v*h/2
            end
        end
        u[:, n+1]=zg(A, H1)
    end
    for n=r-1:N
        gp = [n r a]
        g1 = ones(1, n)
        g2 = g1*wj(gp)
        g=zeros(n,1)
        for i8=1:n
            g[i8, 1]=g2[1, i8]
        end
        u1=zeros(i-1, n)
        for i9=1:n
            u1[:,i9] = u[:, i9]
        end
        H1=-miu*(u1*g)+h^2*f[:, n] + H[:, n]

        A=zeros(i-1,i-1)
        for i10=1:i-1
            A[i10, i10]=miu*g2[1, n+1] + 2*k
            if i10 !== 1
                A[i10, i10-1] = -(k+v*h/2)
            end
            if i10<i-1
                A[i10, i10+1]=-k+v*h/2
            end
        end
        u[:, n+1] = zg(A, H1)
    end
    return u
end

function wj(p)
    n=Int64(p[1, 1]);r=Int64(p[1, 2]);a=p[1, 3]
    A=zeros(n, n+1)
    for iw=1:r-2
        for jw=1:iw+1
        A[iw, jw] = w([iw+1-jw iw+1 iw n a])
        end
    end
    for iw=r-1:n
        for jw=iw-r+2:iw+1
        A[iw, jw]=w([iw+1-jw r iw n a])
        end
    end
    return A
end

function zg(A, H)
    i=length(H)
    r=zeros(i, 1)
    l=zeros(i, 1)
    y=zeros(i, 1)
    U=zeros(i, 1)
    for t=1:i
        if t==1
            r[t, 1] = A[1, 1]
        else
            l[t, 1] = A[t, t-1]/r[t-1, 1]
            r[t, 1] = A[t, t]-l[t, 1]*A[t-1, t]
        end
    end

    for p=1:i
        if p == 1
            y[p, 1] = H[1, 1]
        else
            y[p, 1] = H[p, 1] - l[p, 1]*y[p-1, 1]
        end
    end
    for k=0:i-1
       if k==0
        U[i, 1] = y[i, 1]/r[i, 1]
       else
        U[i-k, 1]=(y[i-k, 1] - A[i-k, i-k+1]*U[i-k+1, 1])/r[i-k, 1]
       end
    end
    return U
end

function w(q)
    i=Int64(q[1, 1]);r=Int64(q[1, 2]);j=q[1, 3];n=q[1, 4];a=q[1, 5]

    ar=ones(1, r-1)
    br=ones(1, r-1)
    for lj=1:r-2
    jj=r-lj-1
    kj=r-2
    tj=i-1
    aj=collect(Int64, 0:tj)
    bj=collect(Int64, tj+2:kj+1)
    cj=[aj; bj]
    dj=collect(Int64, -1:tj-1)
    ej=collect(Int64, tj+1:kj)
    fj=[dj; ej]
    yj = binomial.(cj, jj)
    pj = binomial.(fj, jj)
    sj=1
    tj=1
    for m=1:jj
        sj = sj.*yj[:, m]
        tj = tj.*pj[:, m]
        ar[1, lj] = sum(sj)
        br[1, lj] = sum(tj)
    end
    end
    s=0
    for l=1:r-1
        k=1
        for m1=1:l
            k=k/(m1-a)
        end
        k = k*factorial(l)
        k = k*(ar[1, l]*((n-j)^(l-a))-br[1, l]*((n-j+1)^(l-a)))
        s = s+k
    end
    s = s*((-1)^(i+1))/(factorial(i)*factorial(r-1-i))
    return s
end