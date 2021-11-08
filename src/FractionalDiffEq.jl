module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices, MittagLeffler

include("main.jl")
include("matrix.jl")

export FractionalDiffEqAlgorithm
export PECE, MatrixDiscrete

# Export some special equtions
export bagleytorvik

export solve, FDEProblem

end
