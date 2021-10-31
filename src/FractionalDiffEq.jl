module FractionalDiffEq

using SpecialFunctions, InvertedIndices, MittagLeffler

include("main.jl")
include("matrix.jl")

export FractionalDiffEqAlgorithm
export PECE, MatrixDiscrete

export solve, FDEProblem

end
