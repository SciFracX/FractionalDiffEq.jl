using SciMLBase, DiffEqBase
abstract type AbstractFDEProblem <: SciMLBase.AbstractDEProblem end
abstract type AbstractFODEProblem{uType, tType, oType, isinplace} <: AbstractFDEProblem end
abstract type AbstractFDDEProblem{uType, tType, oType, lType, isinplace} <: AbstractFDEProblem end
abstract type FDEProblem end
struct MultiTermsFODEProblem <: FDEProblem
    parameters
    orders
    rightfun
    rparameters
    rorders
    u0
    tspan
end

#=MultiTermsFODEProblem constructor, dispatch for closedform problem=#
MultiTermsFODEProblem(parameters, orders, rightfun, u0, T) = MultiTermsFODEProblem(parameters, orders, rightfun, nothing, nothing, u0, T)
#MultiTermsFODEProblem(parameters, orders, rightfun, u0, t0, T) = MultiTermsFODEProblem(parameters, orders, rightfun, nothing, nothing, u0, T)


"""

    SingleTermFODEProblem(f, α, u0, tspan)

Define a single term fractional ordinary differential equation, there are only one fractional differential operator in this problem.
"""
#=
abstract type AbstractTestProblem{uType, tType, isinplace} <: SciMLBase.AbstractODEProblem{uType, tType, isinplace} end
=#


SciMLBase.isinplace(prob::AbstractFODEProblem{uType, tType, oType, iip}) where {uType, tType, oType, iip} = iip
SciMLBase.isinplace(prob::AbstractFDDEProblem{uType, tType, oType, lType, isinplace}) where {uType, tType, oType, lType, isinplace} = iip

struct StandardFODEProblem end

struct FODEProblem{uType, tType, oType, isinplace, P, F, bF, PT, K} <:
        AbstractFODEProblem{uType, tType, oType, isinplace}
    f::F
    order::oType
    u0::uType
    tspan::tType
    p::P
    problem_type::PT
    kwargs::K
    SciMLBase.@add_kwonly function FODEProblem{iip}(f::SciMLBase.AbstractODEFunction, order, u0, tspan,
        p = SciMLBase.NullParameters(),
        problem_type = StandardFODEProblem();
        kwargs...) where {iip}
        _tspan = SciMLBase.promote_tspan(tspan)
        #warn_paramtype(p)
        new{typeof(u0), typeof(_tspan), typeof(order), iip, typeof(p),
            typeof(f), typeof(order),
            typeof(problem_type), typeof(kwargs)}(f, order, u0, _tspan, p,
            problem_type, kwargs)
    end

    function FODEProblem{iip}(f, order, u0, tspan, p = SciMLBase.NullParameters(); kwargs...) where {iip}
        FODEProblem(ODEFunction{iip}(f), order, u0, tspan, p; kwargs...)
    end
end

TruncatedStacktraces.@truncate_stacktrace SingleTermFODEProblem 3 1 2

function FODEProblem(f::SciMLBase.AbstractODEFunction, order, u0, tspan, args...; kwargs...)
    FODEProblem{SciMLBase.isinplace(f, 4)}(f, order, u0, tspan, args...; kwargs...)
end

function FODEProblem(f, order, u0, tspan, p = SciMLBase.NullParameters(); kwargs...)
    FODEProblem(ODEFunction(f), order, u0, tspan, p; kwargs...)
end


"""
    FDDEProblem(f, ϕ, α, τ, tspan)

- `f`: The function describing fractional delay differential equations.
- `ϕ`: History function
Construct a fractional delayed differential equation problem.
"""
struct FDDEProblem <: FDEProblem
    f::Function
    ϕ::Union{Number, Function}
    α::Union{Number, Function}
    τ::Union{Number, AbstractArray, Function}
    tspan::Union{Number, Tuple}
end
#=
#=FDDEProblem constructor=#
function FDDEProblem(f::Function,
                     ϕ::Union{Number, Function},
                     α::Union{Number, Function},
                     τ::Union{Number, AbstractArray, Function},
                     T::Union{Tuple, Number})
    return FDDEProblem(f, ϕ, α, τ, T)
end
=#
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
function FODESystem(f::Function,
                    α::AbstractArray,
                    u0::AbstractArray,
                    tspan::Union{Tuple, Number})
    FODESystem(f, α, u0, tspan, nothing)
end


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
function FFPODEProblem(f::Function,
                       order::Union{AbstractArray, Function},
                       u0::Union{AbstractArray, Number},
                       tspan::Union{Tuple, Number})
    FFPODEProblem(f, order, u0, tspan, nothing)
end

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
function FFEODEProblem(f::Function,
                       order::Union{AbstractArray, Function},
                       u0::Union{AbstractArray, Number},
                       tspan::Union{Tuple, Number})
    FFEODEProblem(f, order, u0, tspan, nothing)
end

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
function FFMODEProblem(f::Function,
                       order::Union{AbstractArray, Function},
                       u0::Union{AbstractArray, Number},
                       tspan::Union{Tuple, Number})
    FFMODEProblem(f, order, u0, tspan, nothing)
end

struct FFMODESystem <: FDEProblem
    f::Function
    order::Union{AbstractArray, Function}
    u0::Union{AbstractArray, Number}
    tspan::Union{Tuple, Number}
    p::Union{AbstractArray, Number, Nothing}
end

# If the there are no parameters, we do this:
function FFMODESystem(f::Function,
                      order::Union{AbstractArray, Function},
                      u0::Union{AbstractArray, Number},
                      tspan::Union{Tuple, Number})
    FFMODESystem(f, order, u0, tspan, nothing)
end

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
struct FractionalDiscreteProblem <: FDEProblem
    fun::Function
    α
    u0
    tspan::Union{Tuple, Real, Nothing}
    p::Union{AbstractArray, Nothing}
end

FractionalDiscreteProblem(fun, α, u0, tspan) = FractionalDiscreteProblem(fun, α, u0, tspan, nothing)
FractionalDiscreteProblem(fun, α, u0) = FractionalDiscreteProblem(fun, α, u0, nothing)

"""
    FractionalDifferenceSystem(f, α, u0)

Define Fractional Difference System with the general constructure:

```math
{^G\\nabla_k^\\alpha x(k+1)}=f(x(k))
```

With given initial condition ``x(i)``.
"""
struct FractionalDiscreteSystem <: FDEProblem
    fun::Function
    α
    u0
    tspan::Union{Tuple, Real, Nothing}
    p::Union{AbstractArray, Number, Nothing}
end

FractionalDiscreteSystem(fun, α, u0, tspan) = FractionalDiscreteSystem(fun, α, u0, tspan, nothing)
FractionalDiscreteSystem(fun, α, u0) = FractionalDiscreteSystem(fun, α, u0, nothing)

################################################################################



