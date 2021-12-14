import FractionalDiffEq.FractionalDiffEqAlgorithm

# Algorithms extending

"""
Closed Form solution algorithm.

## Reference:

Dingyu Xue, Northeastern University, China
ISBN:9787030543981
"""
struct ClosedForm <: FractionalDiffEqAlgorithm end

"""
    solve(parameters, order, lparameters, lorders, u, t)

Use Closed-Form solution to obtain numerical solution at zero initial condition.
"""
function solve(parameters, orders, rparameters, rorders, u, t, ::ClosedForm)
    h = t[2]-t[1]
    D = sum(parameters./(h.^orders))
    nT = length(t)

    #TODO: For small orders? e.g. if rparameters and rorders is Number???
    D1 = rparameters[:]./h.^rorders[:]

    nA = length(parameters)
    vect = hcat(orders, rorders)

    y1 = zeros(nT)
    W = ones(nT, length(vect))
    
    @fastmath @inbounds @simd for j=2:nT
        W[j, :] = W[j-1, :].*(1 .-(vect'.+1)/(j.-1))
    end

    @fastmath @inbounds @simd for i=2:nT
        A=y1[i-1:-1:1]'*W[2:i, 1:nA]
        y1[i] = (u[i]-sum(A.*parameters./(h.^orders)))/D
    end

    y = zeros(nT) # Initial the final result

    @fastmath @inbounds @simd for i=2:nT
        y[i] = (W[1:i, (nA+1):end]*D1)'*y1[i:-1:1]
    end

    return y
end