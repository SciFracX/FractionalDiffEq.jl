# Algorithms extending

"""
# Usage

    solve(prob::MultiTermsFODEProblem, ClosedForm())

Use Closed-Form solution to obtain numerical solution at zero initial condition.

## Reference:

Dingyu Xue, Northeastern University, China
ISBN:9787030543981
"""
struct ClosedForm <: FractionalDiffEqAlgorithm end


function solve(prob::MultiTermsFODEProblem, h, ::ClosedForm)
    @unpack parameters, orders, rightfun, rparameters, rorders, u0, tspan = prob
    t0 = tspan[1]; T = tspan[2]
    D = sum(parameters./(h.^orders))
    t = collect(t0:h:T)
    nT = length(t)
    u = rightfun.(t)

    #TODO: For small orders? e.g. if rparameters and rorders is Number???
    isa(rparameters, Number) && typeof(rorders) <: Number ? D1 = rparameters/h^rorders : D1 = rparameters[:]./h.^rorders[:]

    nA = length(parameters)
    vect = hcat(orders, rorders)

    y1 = zeros(nT)
    W = ones(nT, length(vect))
    
    @fastmath @inbounds @simd for j=2:nT
        W[j, :] = W[j-1, :].*(1 .-(vect'.+1)/(j.-1))
    end

    @fastmath @inbounds @simd for i=2:nT
        A = y1[i-1:-1:1]'*W[2:i, 1:nA]
        y1[i] = (u[i]-sum(A.*parameters./(h.^orders)))/D
    end

    y = zeros(nT) # Prelocate the final result

    @fastmath @inbounds @simd for i=2:nT
        y[i] = (W[1:i, (nA+1):end]*D1)'*y1[i:-1:1]
    end
    return FODESolution(t, y)
end