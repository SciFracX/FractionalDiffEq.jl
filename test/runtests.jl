using FractionalDiffEq
using Test

@testset "FractionalDiffEq.jl" begin
    include("test.jl")
    include("auxillary.jl")
    include("models.jl")
end
