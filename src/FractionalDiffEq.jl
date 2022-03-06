module FractionalDiffEq

using LinearAlgebra, SpecialFunctions
using ApproxFun
using InvertedIndices
using QuadGK
using SpecialMatrices
using Polynomials
using FFTW
using UnPack

include("PECE.jl")
include("matrix.jl")
include("MTrap.jl")
include("FLMM.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

include("GL/GLWithMemory.jl")
include("Direct.jl")

include("NonLinear/NonLinear.jl")
include("NonLinear/LorenzADM.jl")

# Fractional partial differential equations
include("FPDE/FiniteDiffEx.jl")
include("FPDE/FiniteDiffIm.jl")
include("FPDE/CaputoDiscrete.jl")

# Fractional delay differential equations
include("FDDE/DelayPECE.jl")
include("FDDE/DelayPI.jl")
include("FDDE/Matrix.jl")
include("FDDE/DelayABM.jl")

# Distributed order differential equations
include("DistributedOrder/utils.jl")
include("DistributedOrder/matrix.jl")

# Fractional Differences equations
include("FractionalDifferences/PECE.jl")

# Fractional integral equations
include("FIE/Qmat.jl")
include("FIE/mycoeffs.jl")

# Mittag Leffler function
include("mlfun.jl")


export FractionalDiffEqAlgorithm

# Export general types
export solve, FDEProblem, FPDEProblem, FDDEProblem

# Detailed problem types
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem, FractionalDifferenceProblem

# FODE solvers
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, G2Direct

# FPDE solvers
export FPDEMatrixDiscrete, FiniteDiffEx, FiniteDiffIm, ADV_DIF

# System of FDE solvers
export NonLinearAlg, GLWithMemory, ModifiedTrap, LorenzADM

# FDDE solvers
export DelayPECE, DelayPI, MatrixForm, DelayABM

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
export PECEDifference




# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special models
export bagleytorvik, diffusion, FractionalLorenz

# Auxiliary functions
export mittleff, mittlefferr, mittleffderiv

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

# Export auxillary functions for FLMM
export ourfft, ourifft, Weights, FastConv

end