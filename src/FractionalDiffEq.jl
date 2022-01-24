module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices

include("PECE.jl")
include("matrix.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

#include("GLWithMemory.jl")
include("Direct.jl")

include("NonLinear/NonLinear.jl")
include("NonLinear/DelayPECE.jl")
include("NonLinear/Qi.jl")

include("mlfun.jl")




export mittleff

export FractionalDiffEqAlgorithm

# Export detailed problem types
export FDEProblem, FPDEProblem, FDDEProblem

export SingleTermFODEProblem, MultiTermsFODEProblem

export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete, ClosedForm, ClosedFormHankelM, G2Direct

export DelayPECE

# Export some api to construct the equation
export RieszMatrix, omega, meshgrid

# Export some special equtions
export bagleytorvik
export diffusion

export NonLinearAlg

export solve, FDEProblem, FODEProblem, FPDEProblem

end