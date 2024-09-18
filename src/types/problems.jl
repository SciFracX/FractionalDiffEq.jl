using SciMLBase
abstract type AbstractFDEProblem <: SciMLBase.AbstractDEProblem end
abstract type FDEProblem end

"""
Defines an multiple terms linear fractional ordinary differential equation (FODE) problem.

## Mathematical Specification of an multi-terms FODE problem

To define an multi-terms FODE Problem, you simply need to given the parameters, their correspoding orders, right hand side function and the initial condition ``u_0``  which define an FODE:

```math
\\frac{du^\\alpha}{d^{\\alpha}t} = f(u,p,t)
```

Multiple terms fractional order differential equations.
"""
struct MultiTermsFODEProblem{uType, tType, oType, pType, F, P, K, isinplace} <:
       SciMLBase.AbstractODEProblem{uType, tType, isinplace}
    parameters::pType
    orders::oType
    f::F
    rparameters::Union{Nothing, pType}
    rorders::Union{Nothing, oType}
    u0::uType
    tspan::tType
    p::P
    kwargs::K

    SciMLBase.@add_kwonly function MultiTermsFODEProblem{iip}(
            parameters, orders, f::SciMLBase.AbstractODEFunction, rparameters,
            rorders, u0, tspan, p = SciMLBase.NullParameters(); kwargs...) where {iip}
        _tspan = SciMLBase.promote_tspan(tspan)
        new{typeof(u0), typeof(_tspan), typeof(orders),
            typeof(parameters), typeof(f), typeof(p), typeof(kwargs), iip}(
            parameters, orders, f, rparameters, rorders, u0, _tspan, p, kwargs)
    end

    function MultiTermsFODEProblem{iip}(
            parameters, orders, f, rparameters, rorders, u0, tspan,
            p = SciMLBase.NullParameters(); kwargs...) where {iip}
        MultiTermsFODEProblem(parameters, orders, ODEFunction{iip}(f),
            rparameters, rorders, u0, tspan, p; kwargs...)
    end
end

function MultiTermsFODEProblem(
        parameters, orders, f, u0, tspan, p = SciMLBase.NullParameters(); kwargs...)
    return MultiTermsFODEProblem{false}(
        parameters, orders, ODEFunction(f), nothing, nothing, u0, tspan, p; kwargs...)
end

"""
Defines an distributed order fractional ordinary differential equation (DODE) problem.

## Mathematical Specification of an distributed order FODE problem

To define an multi-terms FODE Problem, you simply need to given the parameters, their correspoding orders, right hand side function and the initial condition ``u_0``  which define an FODE:

```math
\\frac{du^\\alpha}{d^{\\alpha}t} = f(u,p,t)
```

Distributed order differential equations.
"""
struct DODEProblem{uType, tType, oType, pType, F, P, K, isinplace} <:
       SciMLBase.AbstractODEProblem{uType, tType, isinplace}
    parameters::pType
    orders::oType
    f::F
    u0::uType
    tspan::tType
    p::P
    kwargs::K

    SciMLBase.@add_kwonly function DODEProblem{iip}(
            parameters, orders, f::SciMLBase.AbstractODEFunction, rparameters,
            rorders, u0, tspan, p = SciMLBase.NullParameters(); kwargs...) where {iip}
        _tspan = SciMLBase.promote_tspan(tspan)
        new{typeof(u0), typeof(_tspan), typeof(orders),
            typeof(parameters), typeof(f), typeof(p), typeof(kwargs), iip}(
            parameters, orders, f, rparameters, rorders, u0, _tspan, p, kwargs)
    end

    function DODEProblem{iip}(parameters, orders, f, u0, tspan,
            p = SciMLBase.NullParameters(); kwargs...) where {iip}
        DODEProblem(parameters, orders, ODEFunction{iip}(f), u0, tspan, p; kwargs...)
    end
end

function DODEProblem(
        parameters, orders, f, u0, tspan, p = SciMLBase.NullParameters(); kwargs...)
    return DODEProblem{false}(parameters, orders, ODEFunction(f), u0, tspan, p; kwargs...)
end

struct StandardFODEProblem end

"""
Defines an fractional ordinary differential equation (FODE) problem.

## Mathematical Specification of an FODE problem

To define an FODE Problem, you simply need to given the function ``f`` and the initial condition ``u_0``  which define an FODE:

```math
\\frac{du^\\alpha}{d^{\\alpha}t} = f(u,p,t)
```

There are two different ways of specifying `f`:

  - `f(du,u,p,t)`: in-place. Memory-efficient when avoiding allocations. Best option for most cases unless mutation is not allowed.
  - `f(u,p,t)`: returning `du`. Less memory-efficient way, particularly suitable when mutation is not allowed (e.g. with certain automatic differentiation packages such as Zygote).
  - `order`: the fractional order of the differential equations, commensurate and non-commensurate is both supported.
    -`u₀` should be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
    Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
    provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

`FODEProblem` can be constructed by first building an `ODEFunction` or
by simply passing the FODE right-hand side to the constructor. The constructors
are:

  - `FODEProblem(f::ODEFunction,u0,tspan,p=NullParameters();kwargs...)`
  - `FODEProblem{isinplace,specialize}(f,u0,tspan,p=NullParameters();kwargs...)` :
    Defines the FODE with the specified functions. `isinplace` optionally sets whether
    the function is inplace or not. This is determined automatically, but not inferred.
    `specialize` optionally controls the specialization level. See the
    [specialization levels section of the SciMLBase documentation](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#Specialization-Levels)
    for more details. The default is `AutoSpecialize`.

For more details on the in-place and specialization controls, see the ODEFunction
documentation.

Parameters are optional, and if not given, then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the `ODEFunction` documentation.

### Fields

  - `f`: The function in the ODE.
  - `order`: The order of the FODE.
  - `u0`: The initial condition.
  - `tspan`: The timespan for the problem.
  - `p`: The parameters.
  - `kwargs`: The keyword arguments passed onto the solves.

## Example Problem

```julia
using SciMLBase
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
order = [0.96; 0.96; 0.96]
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = FODEProblem(lorenz!, u0, tspan)

# Test that it worked
using FractionalDiffEq
sol = solve(prob, PIEX())
using Plots;
plot(sol, vars = (1, 2, 3));
```
"""
struct FODEProblem{uType, tType, oType, isinplace, P, F, bF, PT, K} <:
       SciMLBase.AbstractODEProblem{uType, tType, isinplace}
    f::F
    order::oType
    u0::uType
    tspan::tType
    p::P
    problem_type::PT
    kwargs::K
    SciMLBase.@add_kwonly function FODEProblem{iip}(f::SciMLBase.AbstractODEFunction, order,
            u0, tspan, p = SciMLBase.NullParameters(),
            problem_type = StandardFODEProblem(); kwargs...) where {iip}
        _tspan = SciMLBase.promote_tspan(tspan)
        #warn_paramtype(p)
        new{typeof(u0), typeof(_tspan), typeof(order), iip, typeof(p),
            typeof(f), typeof(order), typeof(problem_type), typeof(kwargs)}(
            f, order, u0, _tspan, p, problem_type, kwargs)
    end

    function FODEProblem{iip}(
            f, order, u0, tspan, p = SciMLBase.NullParameters(); kwargs...) where {iip}
        FODEProblem(ODEFunction{iip}(f), order, u0, tspan, p; kwargs...)
    end
end

function FODEProblem(f::SciMLBase.AbstractODEFunction, order, u0, tspan, args...; kwargs...)
    FODEProblem{SciMLBase.isinplace(f, 4)}(f, order, u0, tspan, args...; kwargs...)
end

function FODEProblem(f, order, u0, tspan, p = SciMLBase.NullParameters(); kwargs...)
    FODEProblem(ODEFunction(f), order, u0, tspan, p; kwargs...)
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

struct StandardFDDEProblem end

"""
Defines a fractional delay differential equation (FDDE) problem.

## Mathematical Specification of an FDDE problem

To define an FDDE problem, you simply need to given the function ``f``, the history function ``ϕ``, the fractional order ``α``, the delay ``τ`` and the initial condition ``u₀``  which define an FDDE:

```math
D^\\alpha_t y(t)=f(t,y(t),y(t-τ)),t ≥ 0
```

```math
y(t)=ϕ(t), t ≤ 0
```

There are two different ways of specifying `f`:

  - `f(du,u,h,p,t)`: The function describing fractional delay differential equations.
  - `f(u,h,p,t)`: returning `du`. Less memory-efficient way, particularly suitable when mutation is not allowed (e.g. with certain automatic differentiation packages such as Zygote).
  - `ϕ`: History function
  - `α`: The fractional order of the differential equations, commensurate and non-commensurate is both supported.
  - `τ`: The time delay of the differential equations.

## Problem Type

### Constructors

`FODEProblem` can be constructed by first building an `ODEFunction` or
by simply passing the FODE right-hand side to the constructor. The constructors
are:

  - `FDDEProblem(f::ODEFunction,order,u0,ϕ,tspan,p=NullParameters();kwargs...)`
  - `FDDEProblem{isinplace,specialize}(f,order,u0,ϕ,tspan,p=NullParameters();kwargs...)` :
    Defines the FODE with the specified functions. `isinplace` optionally sets whether
    the function is inplace or not. This is determined automatically, but not inferred.
    `specialize` optionally controls the specialization level. See the
    [specialization levels section of the SciMLBase documentation](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#Specialization-Levels)
    for more details. The default is `AutoSpecialize`.

For more details on the in-place and specialization controls, see the ODEFunction
documentation.

Parameters are optional, and if not given, then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the `ODEFunction` documentation.

### Fields

  - `f`: The function in the FDDE.
  - `ϕ`: The history function in the FDDE.
  - `order`: The order of the FDDE.
  - `τ`: The time delay of the FDDE.
  - `tspan`: The timespan for the problem.
  - `p`: The parameters.
  - `kwargs`: The keyword arguments passed onto the solves.

## Example Problem

```julia
function ϕ(p, t)
    if t == 0
        return 19.00001
    else
        return 19.0
    end
end
function f(y, ϕ, p, t)
    return 3.5 * y * (1 - ϕ / 19)
end
τ = [0.8]
order = 0.97
u0 = 19.00001
tspan = (0.0, 2.0)
dt = 0.5
prob = FDDEProblem(f, order, u0, ϕ, constant_lags = τ, tspan)
sol = solve(prob, DelayPIEX(), dt = dt)
```

"""
struct FDDEProblem{uType, tType, oType, lType, isinplace, P, F, H, K, PT} <:
       SciMLBase.AbstractDDEProblem{uType, tType, lType, isinplace}
    f::F
    order::oType
    u0::uType
    h::H
    tspan::tType
    p::P
    constant_lags::lType
    kwargs::K
    problem_type::PT

    SciMLBase.@add_kwonly function FDDEProblem{iip}(
            f::SciMLBase.AbstractDDEFunction{iip}, order, u0, h,
            tspan, p = SciMLBase.NullParameters(); constant_lags = (),
            problem_type = StandardFDDEProblem(), kwargs...) where {iip}
        _tspan = SciMLBase.promote_tspan(tspan)
        SciMLBase.warn_paramtype(p)
        new{typeof(u0), typeof(_tspan), typeof(order), typeof(constant_lags), isinplace(f),
            typeof(p), typeof(f), typeof(h), typeof(kwargs), typeof(problem_type)}(
            f, order, u0, h, _tspan, p, constant_lags, kwargs, problem_type)
    end

    function FDDEProblem{iip}(f::SciMLBase.AbstractDDEFunction{iip}, h, tspan::Tuple,
            p = SciMLBase.NullParameters(); kwargs...) where {iip}
        FDDEProblem{iip}(f, h(p, first(tspan)), h, tspan, p; kwargs...)
    end
end

FDDEProblem(f, args...; kwargs...) = FDDEProblem(DDEFunction(f), args...; kwargs...)

function FDDEProblem(f::SciMLBase.AbstractDDEFunction, args...; kwargs...)
    FDDEProblem{SciMLBase.isinplace(f)}(f, args...; kwargs...)
end

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
function FODESystem(
        f::Function, α::AbstractArray, u0::AbstractArray, tspan::Union{Tuple, Number})
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
function FFPODEProblem(f::Function, order::Union{AbstractArray, Function},
        u0::Union{AbstractArray, Number}, tspan::Union{Tuple, Number})
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
function FFEODEProblem(f::Function, order::Union{AbstractArray, Function},
        u0::Union{AbstractArray, Number}, tspan::Union{Tuple, Number})
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
function FFMODEProblem(f::Function, order::Union{AbstractArray, Function},
        u0::Union{AbstractArray, Number}, tspan::Union{Tuple, Number})
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
function FFMODESystem(f::Function, order::Union{AbstractArray, Function},
        u0::Union{AbstractArray, Number}, tspan::Union{Tuple, Number})
    FFMODESystem(f, order, u0, tspan, nothing)
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

function FractionalDiscreteProblem(fun, α, u0, tspan)
    FractionalDiscreteProblem(fun, α, u0, tspan, nothing)
end
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

function FractionalDiscreteSystem(fun, α, u0, tspan)
    FractionalDiscreteSystem(fun, α, u0, tspan, nothing)
end
FractionalDiscreteSystem(fun, α, u0) = FractionalDiscreteSystem(fun, α, u0, nothing)

################################################################################
