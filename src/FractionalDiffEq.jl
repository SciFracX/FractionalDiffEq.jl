module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices, FractionalCalculus

include("main.jl")
include("matrix.jl")

export FractionalDiffEqAlgorithm
export PECE, MatrixDiscrete

# Export some api to construct the equation
export D

# Export some special equtions
export bagleytorvik

export solve, FODEProblem, FPDEProblem

end
