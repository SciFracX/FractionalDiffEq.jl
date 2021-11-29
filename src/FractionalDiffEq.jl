module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices

include("main.jl")
include("matrix.jl")
include("Closedform.jl")

export FractionalDiffEqAlgorithm
export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete, ClosedForm

# Export some api to construct the equation
export D, RieszMatrix, omega

# Export some special equtions
export bagleytorvik
export diffusion

export solve, FDEProblem, FODEProblem, FPDEProblem

end