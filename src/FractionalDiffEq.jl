module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices, FractionalCalculus

include("main.jl")
include("matrix.jl")

export FractionalDiffEqAlgorithm
export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete

# Export some api to construct the equation
export D

# Export some special equtions
export bagleytorvik
export diffusion

export solve, FODEProblem, FPDEProblem

end
