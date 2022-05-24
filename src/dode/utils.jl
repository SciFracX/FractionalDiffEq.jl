function DOB(ϕ, alpharange, alphastep, tN, tstep)
    alphas = collect(alpharange[1]:alphastep:alpharange[2])
    alphacount = length(alphas)

    result = zeros(tN, tN)
    phi = ϕ.(alphas)

    for k=1:alphacount
        result = result .+ phi[k]*alphastep*D(tN, alphas[k], tstep)
    end
    return result
end

function DOF(ϕ, alpharange, alphastep, tN, tstep)
    alphas = collect(alpharange[1]:alphastep:alpharange[2])
    alphacount = length(alphas)

    result = zeros(tN, tN)

    phi = ϕ.(alphas)

    for k=1:alphacount
        result = result .+ phi[k]*alphastep*F(tN, alphas[k], tstep)
    end
    return result
end

function DORANORT(ϕ, alpharange, alphastep, tN, tstep)
    alphas = collect(alpharange[1]:alphastep:alpharange[2])
    alphacount = length(alphas)

    result = zeros(tN, tN)

    phi = ϕ.(alphas)

    for k=1:alphacount
        result = result .+ phi[k]*alphastep*ranort(alphas[k], tN, tstep)
    end
    return result
end

function ranort(alpha, N)
    k=collect(0:N-1)
    rc = ((-1)*ones(size(k))).^k.*gamma.(alpha+1).*(gamma.(alpha*0.5 .-k.+1).*gamma.(alpha*0.5 .+ k.+1)).^(-1)
    rc = rc*(cos(alpha*π*0.5))

    R = zeros(N, N)

    for m=1:N
        R[m, m:N] = rc[1:N-m+1]
    end

    for i=1:N-1
        for j=i:N
            R[j, i] = R[i, j]
        end
    end
    return R
end
# Multiple dispatch for ranort
function ranort(alpha, N, h)
    R = ranort(alpha, N)
    R = R*h^(-alpha)
    return R
end