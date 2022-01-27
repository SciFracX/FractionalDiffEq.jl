import FractionalDiffEq: FDEProblem, SingleTermFODEProblem, MultiTermsFODEProblem

using SpecialFunctions


"""
    FDEProblem

General parent type for all kinds of problems in FractionalDiffEq.jl.
"""
abstract type FDEProblem end

"""
    MultiTermsFODEProblem(parameters, orders, rparameters, rorders)

Define a multi-terms fractional ordinary differential equation.
"""
struct MultiTermsFODEProblem <: FDEProblem
    parameters
    orders
    rightfun
    rparameters::Union{Matrix, Nothing}
    rorders::Union{Matrix, Nothing}
end

#=MultiTermsFODEProblem constructor=#
MultiTermsFODEProblem(parameters, orders, rightfun) = MultiTermsFODEProblem(parameters, orders, rightfun, nothing, nothing)


"""

    SingleTermFODEProblem(f, α, h)

Define a single term fractional ordinary differential equation, there are only one term in this problem.
"""
struct SingleTermFODEProblem <: FDEProblem
    f
    α
    u0
    T
end


"""
Base type of FractionalDiffEq algorithms
"""
abstract type FractionalDiffEqAlgorithm end


"""
Predict-Evaluate-Correct-Evaluate algorithm.

For more details, please refer to [Predictor-Corrector algorithms](https://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method)

This PECE algorithm is taken from Diethelm's paper.

```tex
@article{
title={A predictor-corrector approach for the numerical solution of fractional differential equations},
author={Diethelm, Kai and Ford, Neville J. and Freed, Alan D.}
doi={https://doi.org/10.1023/A:1016592219341}
}
```
"""
struct PECE <: FractionalDiffEqAlgorithm end
#TODO: Use Richardson extrapolation to refine the PECE algorithms 






"""

    FPDEProblem(α, β, T, M, N)

"""
struct FPDEProblem <: FDEProblem
    α
    β
    T
    M
    N
end

"""
    FDDEProblem(f, ϕ, α, τ)

Construct a fractional delayed differential equation peoblem.
"""
struct FDDEProblem <: FDEProblem
    f
    ϕ
    α
    τ
    t0::Union{Number, Nothing}
end

#=FDDEProblem constructor=#
FDDEProblem(f, ϕ, α, τ) = FDDEProblem(f, ϕ, α, τ, nothing)


"""
    solve(FODEProblem, PECE())

After define the FDEProblem, use **PECE(Predict-Evaluate-Correct-Evaluate) algorithm** to computing the Fractional Differential Equation
"""
function solve(FODE::SingleTermFODEProblem, h, ::PECE)
    f, α, u0, T = FODE.f, FODE.α, FODE.u0, FODE.T
    N = Int64(floor(T/h))
    y = zeros(N+1)
    leftsum = zero(Float64)
    l = floor(α)

    # Handling initial value
    if l == 0
        leftsum = u0
    elseif l == 1
        leftsum = u0 + T*u0
    end

    @fastmath @inbounds @simd for n ∈ 0:N
        y[n+1] = leftsum + h^α/gamma(α+2)*(f((n+1)*h, predictor(f, y, α, n, h, u0, T)) + right(f, y, α, n, h))
    end

    return y
end

function right(f, y, α, n, h)
    temp = zero(Float64)

    @fastmath @inbounds @simd for j = 0:n
        temp += A(j, n, α)*f(j*h, y[Int64(j+1)])
    end

    return temp
end

function predictor(f, y, α::Float64, n::Int64, h, u0, T)
    predict = 0
    leftsum = 0

    l = floor(α)

    # Handling initial value
    if l == 0
        leftsum = u0
    elseif l == 1
        leftsum = u0 + T*u0
    end

    @fastmath @inbounds @simd for j ∈ 0:n
        predict += B(j, n, α, h)*f(j*h, y[j+1])
    end

    return leftsum + predict
end


function A(j, n, α)
    if j == 0
        return n^(α+1) - (n-α)*(n+1)^α
    elseif 1 ≤ j ≤ n
        return (n-j+2)^(α+1) + (n-j)^(α+1) - 2*(n-j+1)^(α+1)
    elseif j == n+1
        return 1
    end
end

function B(j, n, α, h)
    return h^α/α*((n + 1 - j)^α - (n - j)^α)
end