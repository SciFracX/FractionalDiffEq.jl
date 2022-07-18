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
abstract type FDOAlgorithm <: AbstractFDEAlgorithm end
abstract type FDDEAlgorithm <: AbstractFDEAlgorithm end
abstract type FIEAlgorithm <: AbstractFDEAlgorithm end
abstract type FPDEAlgorithm <: AbstractFDEAlgorithm end
abstract type FODESystemAlgorithm <: AbstractFDEAlgorithm end
abstract type MultiTermsFDEAlgorithm <: AbstractFDEAlgorithm end
abstract type SingleTermFDEAlgorithm <: AbstractFDEAlgorithm end




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

Construct a fractional matrix differential equation with delay with general form:

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
    FFPODEProblem(f, α, u0, tspan, p)

Define fractal-fractional differential equations problems with power law kernel.
"""
struct FFPODEProblem <: FDEProblem
    f::Function
    order::Union{AbstractArray, Function}
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
FFPODEProblem(f, order, u0, tspan) = FFPODEProblem(f, order, u0, tspan, nothing)

"""
    FFEODEProblem(f, α, u0, tspan, p)

Define fractal-fractional differential equations problems with expoenntial decay kernel.
"""
struct FFEODEProblem <: FDEProblem
    f::Function
    order::Union{AbstractArray, Function}
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
FFEODEProblem(f, order, u0, tspan) = FFEODEProblem(f, order, u0, tspan, nothing)

"""
    FFMODEProblem(f, α, u0, tspan, p)

Define fractal-fractional differential equations problems with Mittag Leffler kernel.
"""
struct FFMODEProblem <: FDEProblem
    f::Function
    order::Union{AbstractArray, Function}
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
FFMODEProblem(f, order, u0, tspan) = FFMODEProblem(f, order, u0, tspan, nothing)



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
{^G\\nabla_k^\\alpha x(k+1)}=f(x(k))
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



