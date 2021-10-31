module FractionalDiffEq

using SpecialFunctions, InvertedIndices, MittagLeffler

include("main.jl")

export FractionalDiffEqAlgorithm
export PECE, MatrixDiscrete

export solve, FDEProblem

end
