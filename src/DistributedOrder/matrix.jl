import FractionalDiffEq: FractionalDiffEqAlgorithm, solve, eliminator, DOB

using LinearAlgebra

struct DOMatrixDiscrete <: FractionalDiffEqAlgorithm end

function testsolve(ω, t, h, B)
    N = length(t)+1
    M = zeros(N, N)

    M = DOB(ω, [0, 1], 0.01, N-1, h)+B*(zeros(N-1, N-1)+I)
    #F = ω.(t)

    M = eliminator(N-1, 1)*M*eliminator(N-1, 1)'
    F = eliminator(N-1, 1)*F

    Y = M\F

    Y0 = [0; Y]

    return Y0.+1

end
#=
h=0.01
t=collect(0:h:5)
result = testsolve(x->6*x*(1-x), t, h, 0.1)
#plot(t, result)
=#