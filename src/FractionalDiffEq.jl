module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices

include("main.jl")
include("matrix.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/main.jl")


include("mlfun.jl")
export mittleff

export FractionalDiffEqAlgorithm
export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete, ClosedForm, ClosedFormHankelM

# Export some api to construct the equation
export D, RieszMatrix, omega

# Export some special equtions
export bagleytorvik
export diffusion

export solve, FDEProblem, FODEProblem, FPDEProblem

end