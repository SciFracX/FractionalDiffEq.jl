function Base.show(io::IO, sol::FODESolution)
    show(io::IO, typeof(sol.t))
    show(stdout, MIME("text/plain"), sol.t)
    println()
    show(io::IO, typeof(sol.u))
    show(stdout, MIME("text/plain"), sol.u)
end

function Base.show(io::IO, prob::SingleTermFODEProblem)
    printstyled("SingleTermFODEProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("timespan: $(prob.T)")
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::MultiTermsFODEProblem)
    printstyled("MultiTermsFODEProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.orders)", color=:red)
    println()
    println("timespan: $(prob.tspan)")
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FODESystem)
    printstyled("FODESystem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("timespan: $(prob.T)")
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FDDEProblem)
    printstyled("FDDEProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("timespan: $(prob.T)")
    println("history function: $(prob.ϕ)")
end

function Base.show(io::IO, prob::FDDESystem)
    printstyled("FDDESystem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("timespan: $(prob.T)")
    println("history function: $(prob.ϕ)")
end

function Base.show(io::IO, prob::DODEProblem)
    printstyled("DODEProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.orders)", color=:red)
    println()
    println("timespan: $(prob.tspan)")
end

function Base.show(io::IO, prob::FractionalDifferenceProblem)
    printstyled("FractionalDifferenceProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FIEProblem)
    printstyled("FIEProblem", color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.orders)", color=:red)
    println()
    println("timespan: $(prob.tspan)")
end



"""
Fractional differential equation solutions visulization hooks.
"""
@recipe f(sol::FODESolution) = sol.t, sol.u

@recipe f(sol::FDifferenceSolution) = sol.t, sol.u

@recipe f(sol::FIESolution) = sol.t, sol.u

@recipe f(sol::DODESolution) = sol.t, sol.u