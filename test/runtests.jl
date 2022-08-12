using FractionalDiffEq, SpecialFunctions
using Test

@testset "FractionalDiffEq.jl" begin
    include("FODETests.jl")
    include("FDDETests.jl")
    include("FPDETests.jl")
    include("FDifferenceTests.jl")
    include("DODETests.jl")
    include("FFODETests.jl")
    include("auxillary.jl")
    include("FOLEtest.jl")
end
