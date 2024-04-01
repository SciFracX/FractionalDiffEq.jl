using FractionalDiffEq, SpecialFunctions, SciMLBase
using Test

function test_sol(sol::SciMLBase.AbstractODESolution)
    return transpose(mapreduce(permutedims, vcat, sol.u))
end

function test_sol(u::AbstractArray)
    return transpose(mapreduce(permutedims, vcat, u))
end

@testset "FractionalDiffEq.jl" begin
    include("fode.jl")
    include("fdde.jl")
    include("discrete.jl")
    #include("dode.jl")
    include("ffode.jl")
    include("auxillary.jl")
    #include("fole.jl")
end
