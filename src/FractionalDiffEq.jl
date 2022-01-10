module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices

include("PECE.jl")
include("matrix.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")

include("NonLinear/NonLinear.jl")

include("mlfun.jl")




export mittleff

export FractionalDiffEqAlgorithm

# Export problem types
export FDEProblem, FPDEProblem

export SingleTermFODEProblem, MultiTermsFODEProblem

export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete, ClosedForm, ClosedFormHankelM

# Export some api to construct the equation
export D, RieszMatrix, omega

# Export some special equtions
export bagleytorvik
export diffusion

export NonLinearAlg

export solve, FDEProblem, FODEProblem, FPDEProblem

end