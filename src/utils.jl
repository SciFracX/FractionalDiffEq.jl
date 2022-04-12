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
    println("timespan: $(prob.T)")
    println("u0: $(prob.u0)")
end



"""
Fractional ordinary differential equations solutions visulization hooks.
"""
@recipe f(sol::FODESolution) = sol.t, sol.u