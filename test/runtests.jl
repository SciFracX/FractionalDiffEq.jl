using FractionalDiffEq
using Test
using MittagLeffler

@testset "FractionalDiffEq.jl" begin
    include("test.jl")
    include("auxillary.jl")
end
