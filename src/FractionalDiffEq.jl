module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, InvertedIndices, QuadGK, SpecialMatrices

include("PECE.jl")
include("matrix.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

include("GL/GLWithMemory.jl")
include("Direct.jl")

include("NonLinear/NonLinear.jl")
include("NonLinear/Qi.jl")

# Fractional delay differential equations
include("FDDE/DelayPECE.jl")
include("FDDE/PI.jl")

# Auxiliary functions
include("mlfun.jl")




export mittleff, mittlefferr, mittleffderiv

export FractionalDiffEqAlgorithm

# Export detailed problem types
export FDEProblem, FPDEProblem, FDDEProblem

export SingleTermFODEProblem, MultiTermsFODEProblem

export PECE, FODEMatrixDiscrete, FPDEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, G2Direct

export GLWithMemory

export DelayPECE, DelayPI

# Export some api to construct the equation
export RieszMatrix, omega, meshgrid

# Export some special equtions
export bagleytorvik
export diffusion

export NonLinearAlg

export solve, FDEProblem, FODEProblem, FPDEProblem

end