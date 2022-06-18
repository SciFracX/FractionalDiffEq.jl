"""
    FDEProblem

General type for all kinds of problems in FractionalDiffEq.jl.
"""
abstract type FDEProblem end

"""
    MultiTermsFODEProblem(parameters, orders, rightfun, u0, T)

    MultiTermsFODEProblem(parameters, orders, rightfun, rparameters, rorders, u0, T)

Define a multi-terms fractional ordinary differential equation.
"""
struct MultiTermsFODEProblem <: FDEProblem
    parameters::AbstractArray
    orders::AbstractArray
    rightfun::Union{Function, Number}
    rparameters::Union{AbstractArray, Nothing}
    rorders::Union{AbstractArray, Nothing}
    u0::AbstractArray
    tspan::Union{Tuple{Number, Number}, Number}
end

#=MultiTermsFODEProblem constructor, dispatch for closedform problem=#
MultiTermsFODEProblem(parameters, orders, rightfun, u0, T) = MultiTermsFODEProblem(parameters, orders, rightfun, nothing, nothing, u0, T)
#MultiTermsFODEProblem(parameters, orders, rightfun, u0, t0, T) = MultiTermsFODEProblem(parameters, orders, rightfun, nothing, nothing, u0, T)


"""

    SingleTermFODEProblem(f, α, u0, tspan)

Define a single term fractional ordinary differential equation, there are only one fractional differential operator in this problem.
"""
struct SingleTermFODEProblem <: FDEProblem
    f::Function
    α::Float64
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
end


"""
Base type for all of the FractionalDiffEq algorithms
"""
abstract type AbstractFDEAlgorithm end



"""
    FPDEProblem(α, β, T, M, N)

Define fractional order partial differential equations problem
"""
struct FPDEProblem <: FDEProblem
    α
    β
    T
    M
    N
end

"""
    FDDEProblem(f, ϕ, α, τ, tspan)

Construct a fractional delayed differential equation problem.
"""
struct FDDEProblem <: FDEProblem
    f::Function
    ϕ::Union{Number, Function}
    α::Union{Number, Function}
    τ::Union{Number, AbstractArray, Function}
    tspan::Union{Number, Tuple}
end

#=FDDEProblem constructor=#
FDDEProblem(f, ϕ, α, τ, T) = FDDEProblem(f, ϕ, α, τ, T, nothing)

"""
    FDDESystem(f, ϕ, α, τ, T)

Construct system of fractional delay differential equations problem.
"""
struct FDDESystem <: FDEProblem
    f::Function
    ϕ::AbstractArray
    α::Union{AbstractArray, Number}
    τ::Number
    T
end

"""
    FDDEMatrixProblem(α, τ, A, B, f, x0, tspan)

Construct a fractional matrix differential equation with delay. The general type is:

```math
D_{t_0}^\\alpha\\textbf{x}(t)=\\textbf{A}(t)\\textbf{x}(t)+\\textbf{B}(t)\\textbf{x}(t-\\tau)+\\textbf{f}(t)
```
"""
struct FDDEMatrixProblem <: FDEProblem
    α::Float64
    τ::Float64
    A::Matrix
    B::Matrix
    f::Vector
    x0::Function
    tspan::Tuple{Float64, Float64}
end

"""
    FODESystem(f, α, u0, tspan, p)

Define system of fractional differential equations.
"""
struct FODESystem <: FDEProblem
    f::Function
    α::AbstractArray
    u0::AbstractArray
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
FODESystem(f, α, u0, tspan) = FODESystem(f, α, u0, tspan, nothing)


"""
    FFODEProblem(f, α, u0, tspan, p)

Define fractal-fractional differential equations problems.
"""
struct FFODEProblem <: FDEProblem
    f::Function
    order::AbstractArray
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
FFODEProblem(f, order, u0, tspan) = FFODEProblem(f, order, u0, tspan, nothing)




"""
    DODEProblem(parameters, orders, interval, tspan, rightfun)

Define distributed order differential equation problem.
"""
struct DODEProblem <: FDEProblem
    parameters::AbstractArray
    orders::AbstractArray
    interval::Tuple
    rightfun::Function
    u0
    tspan
end

"""
    FractionalDifferenceProblem(f, α, x0)

Define fractional difference equation problems.
"""
struct FractionalDifferenceProblem <: FDEProblem
    fun::Function
    α
    u0
end

"""
    FractionalDifferenceSystem(f, α, u0)

Define Fractional Difference System with the general constructure:

```math
{_{a-m}^G\\nabla_k^\\alpha x(k+1)}=f(x(k))\\
```

With given initial condition ``x(i)``.
"""
struct FractionalDifferenceSystem <: FDEProblem
    fun::Function
    α
    u0
end


"""
    FIEProblem(parameters, orders, rightfun, tspan)

Define fractional integral equation problems.
"""
struct FIEProblem <: FDEProblem
    parameters::AbstractArray
    orders::AbstractArray
    rightfun::Union{Function, Number}
    tspan
end


################################################################################

abstract type AbstractFDESolution end

struct FODESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end


struct FDifferenceSolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FIESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct DODESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end


"""
# Usage

    solve(prob::SingleTermFODEProblem, h, PECE())

Predict-Evaluate-Correct-Evaluate algorithm.

For more details, please refer to [Predictor-Corrector algorithms](https://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method)

This PECE algorithm is taken from Diethelm's paper.

### References

```tex
@article{
title={A predictor-corrector approach for the numerical solution of fractional differential equations},
author={Diethelm, Kai and Ford, Neville J. and Freed, Alan D.}
doi={https://doi.org/10.1023/A:1016592219341}
}
```
"""
struct PECE <: AbstractFDEAlgorithm end
#TODO: Use Richardson extrapolation to refine the PECE algorithms 
#TODO: Rename as ABM

"""
    solve(problem::SingleTermFODEProblem, h, PECE())

Generalproblem solving API for solving FDE problems.
"""
function solve(FODE::SingleTermFODEProblem, h::Float64, ::PECE)
    @unpack f, α, u0, tspan = FODE
    t0 = tspan[1]; T = tspan[2]
    N::Int64 = round(Int, (T-t0)/h)
    y = zeros(N+1) # Initialize the solution array


    @fastmath @inbounds @simd for n ∈ 0:N
        # Handling initial value
        temp=0
        for k=0:ceil(Int, α)-1
            temp += u0[k+1]*(t0+(n+1)*h)^k/factorial(k)
        end
        y[n+1] = temp + h^α/gamma(α+2)*(f(t0+(n+1)*h, predictor(f, y, α, n, h, u0, t0)) + right(f, y, α, n, h))
    end

    tspan = collect(t0:h:T)
    return FODESolution(tspan, y)
end

function right(f, y, α, n, h::Float64)
    temp = zero(Float64)

    @fastmath @inbounds @simd for j = 0:n
        temp += A(j, n, α)*f(j*h, y[j+1])
    end

    return temp
end

function predictor(f::Function, y, α::Float64, n::Integer, h::Float64, u0, t0)
    predict = 0
    leftsum = 0

    l = ceil(Int, α)

    # Handling initial value
    for k=0:l-1
        leftsum += u0[k+1]*(t0+(n+1)*h)^k/factorial(k)
    end

    @turbo for j ∈ 0:n
        predict += B(j, n, α)*f(j*h, y[j+1])
    end

    return leftsum + h^α/α*predict
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

B(j, n, α) = ((n + 1 - j)^α - (n - j)^α) # Moved the h^α/α to the end of predictor: return leftsum + h^α/α*predict