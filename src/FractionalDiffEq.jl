module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices

include("main.jl")
include("matrix.jl")

export FractionalDiffEqAlgorithm
export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete

# Export some api to construct the equation
export D, RieszMatrix, omega

# Export some special equtions
export bagleytorvik
export diffusion

export solve, FODEProblem, FPDEProblem

end