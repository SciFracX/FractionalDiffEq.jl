module FractionalDiffEq

using LinearAlgebra
using SpecialFunctions
using InvertedIndices
using QuadGK
using SpecialMatrices
using Polynomials

include("PECE.jl")
include("matrix.jl")
include("MTrap.jl")
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

# Distributed order differential equations
include("DistributedOrder/utils.jl")
include("DistributedOrder/matrix.jl")

# Mittag Leffler function
include("mlfun.jl")


export FractionalDiffEqAlgorithm

# Export general types
export solve, FDEProblem, FPDEProblem, FDDEProblem

# Detailed problem types
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem

# FODE solvers
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, G2Direct

# FPDE solvers
export FPDEMatrixDiscrete

# System of FDE solvers
export NonLinearAlg, GLWithMemory, ModifiedTrap

# FDDE solvers
export DelayPECE, DelayPI

# DODE solvers
export DOMatrixDiscrete

# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special equtions
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittlefferr, mittleffderiv

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isfunction

end