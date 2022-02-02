import FractionalDiffEq: D, F

function DOB(ϕ, alpharange, alphastep, tN, tstep)
    alphas = collect(alpharange[1]:alphastep:alpharange[2])
    alphacount = length(alphas)

    result = zeros(tN, tN)

    phi = ϕ.(alf)

    for k=1:alphacount
        result = result .+ phi[k]*alphastep*D(tN, alphas[k], tstep)
    end
    return result
end

function DOF(ϕ, alpharange, alphastep, tN, tstep)
    alphas = collect(alpharange[1]:alphastep:alpharange[2])
    alphacount = length(alphas)

    result = zeros(tN, tN)

    phi = ϕ.(alf)

    for k=1:alphacount
        result = result .+ phi[k]*alphastep*F(tN, alphas[k], tstep)
    end
    return result
end