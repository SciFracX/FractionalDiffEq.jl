using FractionalDiffEq
using Test

@testset "FractionalDiffEq.jl" begin
    include("FODETests.jl")
    include("FDDETests.jl")
    include("FPDETests.jl")
    include("FDifferenceTests.jl")
    include("FIETests.jl")
    include("DODETests.jl")
    include("auxillary.jl")
    include("models.jl")
end
