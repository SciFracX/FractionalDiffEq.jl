function Base.show(io::IO, sol::FODESolution)
    println("Time span $(typeof(sol.t))")
    println(sol.t)
    println("Solution $(typeof(sol.u))")
    println(sol.u)
end

function Base.show(io::IO, prob::MultiTermsFODEProblem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.orders)", color=:red)
    println()
    println("timespan: $(prob.tspan)")
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FDDESystem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("timespan: $(prob.T)")
    println("history function: $(prob.ϕ)")
end


function Base.show(io::IO, prob::FractionalDiscreteProblem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.α)", color=:red)
    println()
    println("u0: $(prob.u0)")
end

function Base.show(io::IO, prob::FFPODEProblem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.order[1])", color=:red)
    printstyled(" and $(prob.order[2])", color=:red)
    println()
    println("timespan: $(prob.tspan)")
end

function Base.show(io::IO, prob::FFEODEProblem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.order[1])", color=:red)
    printstyled(" and $(prob.order[2])", color=:red)
    println()
    println("timespan: $(prob.tspan)")
end

function Base.show(io::IO, prob::FFMODEProblem)
    printstyled(typeof(prob), color=:light_blue)
    printstyled(" with order ")
    printstyled("$(prob.order[1])", color=:red)
    printstyled(" and $(prob.order[2])", color=:red)
    println()
    println("timespan: $(prob.tspan)")
end


function Base.show(io::IO, LE::FOLE)
    printstyled("Fractional Lyapunov exponents:", color=:light_blue)
    printstyled("$(LE.LE)")
    println()
    printstyled("Timespan:", color=:light_blue)
    printstyled("$(LE.t)")
end


"""
Fractional differential equation solutions visulization hooks.
"""
@recipe f(sol::FODESolution) = sol.t, sol.u

@recipe f(sol::FDifferenceSolution) = sol.t, sol.u

@recipe f(sol::DODESolution) = sol.t, sol.u

@recipe f(sol::FFMODESolution) = sol.t, sol.u

@recipe function f(sol::FODESystemSolution; vars=nothing)
    if typeof(vars) == Nothing # When vars is not specified, show time versus each variable.
        l = size(sol.u, 1)
        for i in 1:l
            @series begin
                sol.t, sol.u[i, :]
            end
        end
    else
        index = Int[]
        for i in vars
            append!(index, i)
        end
        len = length(index)
        index0 = findall(x->x==0, index) # Find the 0 index, which is the time index.
        if length(index0) == 0
            if len == 3
                sol.u[index[1], :], sol.u[index[2], :], sol.u[index[3], :]
            elseif len == 2
                sol.u[index[1], :], sol.u[index[2], :]
            else
                error("Plots variable index is not correct.")
            end
        elseif length(index0) == 1
            newindex = deleteat!(index, index0)
            newlen = length(newindex)
            if newlen == 2
                for i in newindex
                    @series begin
                        sol.t, sol.u[i, :]
                    end
                end
            elseif newlen == 1
                sol.t, sol.u[newindex[1], :]
            end
        end
    end
end


@recipe function f(sol::FDDESystemSolution; vars=nothing)
    if typeof(vars) == Nothing # When vars is not specified, show time versus each variable.
        l = size(sol.u, 1)
        for i in 1:l
            @series begin
                sol.t, sol.u[i, :]
            end
        end
    else
        index = Int[]
        for i in vars
            append!(index, i)
        end
        len = length(index)
        index0 = findall(x->x==0, index) # Find the 0 index, which is the time index.
        if length(index0) == 0
            if len == 3
                sol.u[index[1], :], sol.u[index[2], :], sol.u[index[3], :]
            elseif len == 2
                sol.u[index[1], :], sol.u[index[2], :]
            else
                error("Plots variable index is not correct.")
            end
        elseif length(index0) == 1
            newindex = deleteat!(index, index0)
            newlen = length(newindex)
            if newlen == 2
                for i in newindex
                    @series begin
                        sol.t, sol.u[i, :]
                    end
                end
            elseif newlen == 1
                sol.t, sol.u[newindex[1], :]
            end
        end
    end
end


@recipe function f(sol::FFMODESystemSolution; vars=nothing)
    if typeof(vars) == Nothing # When vars is not specified, show time versus each variable.
        l = size(sol.u, 1)
        for i in 1:l
            @series begin
                sol.t, sol.u[i, :]
            end
        end
    else
        index = Int[]
        for i in vars
            append!(index, i)
        end
        len = length(index)
        index0 = findall(x->x==0, index) # Find the 0 index, which is the time index.
        if length(index0) == 0
            if len == 3
                sol.u[index[1], :], sol.u[index[2], :], sol.u[index[3], :]
            elseif len == 2
                sol.u[index[1], :], sol.u[index[2], :]
            else
                error("Plots variable index is not correct.")
            end
        elseif length(index0) == 1
            newindex = deleteat!(index, index0)
            newlen = length(newindex)
            if newlen == 2
                for i in newindex
                    @series begin
                        sol.t, sol.u[i, :]
                    end
                end
            elseif newlen == 1
                sol.t, sol.u[newindex[1], :]
            end
        end
    end
end



@recipe function f(LyapunovExponents::FOLE; vars=nothing)
    if typeof(vars) == Nothing # When vars is not specified, show time versus each variable.
        l = size(LyapunovExponents.LE, 1)
        for i in 1:l
            @series begin
                LyapunovExponents.t, LyapunovExponents.LE[i, :]
            end
        end
    else
        index = Int[]
        for i in vars
            append!(index, i)
        end
        len = length(index)
        index0 = findall(x->x==0, index) # Find the 0 index, which is the time index.
        if length(index0) == 0
            if len == 3
                LyapunovExponents.LE[index[1], :], LyapunovExponents.LE[index[2], :], LyapunovExponents.LE[index[3], :]
            elseif len == 2
                LyapunovExponents.LE[index[1], :], LyapunovExponents.LE[index[2], :]
            else
                error("Plots variable index is not correct.")
            end
        elseif length(index0) == 1
            newindex = deleteat!(index, index0)
            newlen = length(newindex)
            if newlen == 2
                for i in newindex
                    @series begin
                        LyapunovExponents.t, LyapunovExponents.LE[i, :]
                    end
                end
            elseif newlen == 1
                LyapunovExponents.t, LyapunovExponents.LE[newindex[1], :]
            end
        end
    end
end