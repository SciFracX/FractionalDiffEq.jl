using FractionalDiffEq, SpecialFunctions
using Test

function test_sol(sol)
    return mapreduce(permutedims, vcat, sol.u)
end

@testset "FractionalDiffEq.jl" begin
    include("fode.jl")
    include("fdde.jl")
    include("discrete.jl")
    include("dode.jl")
    include("ffode.jl")
    include("auxillary.jl")
    include("fole.jl")
end
