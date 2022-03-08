function myeval(u, x, q)
    
    if size(u, 2) > 1
        uu = zeros(length(x), size(u, 2))
        for k = 1:size(u, 2)
            uu[:, k] = myeval(u[:, k], x, q)
        end
        return uu
    end
    
    N::Int64 = length(u)/q
    
    if q == 2
        uu = clenshawP(x, u[1:N]) .+ sqrt.(1 .+x).*clenshawU(x, u[N+1:2*N])
    else
        uu = 0*x
        for k = 0:q-1
            uu = uu .+ (1+x).^(k/q).*clenshawJ_special(x, u[k*N+collect(1:N)], k/q)
        end
    end
    return uu
end
    
function clenshawP(x, c)
    bk1 = 0*x
    bk2 = bk1
    N = size(c, 1)-1
    for k = N:-1:1
        bk = c[k+1] .+ (2*k+1)/(k+1)*x.*bk1 .- (k+1)/(k+2)*bk2
        bk2 = bk1
        bk1 = bk
    end
    y = c[1] .+ x.*bk1 - 0.5*bk2
    return y
end
    
function clenshawU(x, c)
    bk1 = 0*x
    bk2 = bk1
    N = size(c, 1)-1
    for k = N:-1:1
        bk = c[k+1] .+ 2*x.*bk1 - bk2
        bk2 = bk1
        bk1 = bk
    end
    y = c[1] .+ 2*x.*bk1 - bk2
    return y
end
    
function clenshawJ_special(x, c, b)
    bk1 = 0*x
    bk2 = bk1
    N = size(c, 1)-1
    for k = N:-1:1
        Ak = (2*k+3)
        Bk = (1-2*b)/(2*k+1)
        Ck = (k+2-b)*(k+1+b)*(2*k+5)/((k+3)*(2*k+3))
        bk = c[k+1] + ((Ak*x+Bk).*bk1 - Ck*bk2)/(k+2)
        bk2 = bk1
        bk1 = bk
    end
    Ck = (2-b)*(1+b)*5/18
    y = c[1] + 0.5*(3*x+1-2*b).*bk1 - Ck*bk2
    return y
end