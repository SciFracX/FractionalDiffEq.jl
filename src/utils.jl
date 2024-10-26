function Base.show(io::IO, sol::FODESolution)
    println("Time span $(typeof(sol.t))")
    println(sol.t)
    println("Solution $(typeof(sol.u))")
    println(sol.u)
end

function Base.show(io::IO, prob::MultiTermsFODEProblem)
    printstyled(typeof(prob), color = :light_blue)
    printstyled(" with order ")
    printstyled("$(prob.orders)", color = :red)
    println()
    println("timespan: $(prob.tspan)")
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FDDESystem)
    printstyled(typeof(prob), color = :light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color = :red)
    println()
    println("timespan: $(prob.T)")
    println("history function: $(prob.ϕ)")
end

function Base.show(io::IO, prob::FractionalDiscreteProblem)
    printstyled(typeof(prob), color = :light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color = :red)
    println()
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, LE::FOLE)
    printstyled("Fractional Lyapunov exponents:", color = :light_blue)
    printstyled("$(LE.LE)")
    println()
    printstyled("Timespan:", color = :light_blue)
    printstyled("$(LE.t)")
end
