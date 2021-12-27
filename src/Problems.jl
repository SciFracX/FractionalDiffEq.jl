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
    rparameters
    rorders
end


"""

    SingleTermFODEProblem(f, α, h)

Define a single term fractional ordinary differential equation, there are only one term in this problem.
"""
struct SingleTermFODEProblem <: FDEProblem
    f
    α
    h
end