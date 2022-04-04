module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, SparseArrays
using ApproxFun
using InvertedIndices
using QuadGK
using SpecialMatrices
using Polynomials
using FFTW
using UnPack
using LoopVectorization
using HypergeometricFunctions
using ToeplitzMatrices

include("PECE.jl")
include("matrix.jl")
include("PI.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

include("GL/GLWithMemory.jl")
include("GL/GL.jl")
include("IMEX/IMEX_I.jl")
include("IMEX/IMEX_I_A.jl")
include("IMEX/IMEX_II.jl")
#include("Direct.jl")

include("NonLinear/NonLinear.jl")
include("NonLinear/LorenzADM.jl")

# Fractional partial differential equations
include("FPDE/FiniteDiffEx.jl")
include("FPDE/FiniteDiffIm.jl")
include("FPDE/CaputoDiscrete.jl")
include("FPDE/GL.jl")

# Fractional delay differential equations
include("FDDE/DelayPECE.jl")
include("FDDE/DelayPI.jl")
include("FDDE/Matrix.jl")
include("FDDE/DelayABM.jl")
include("FDDE/Yuan.jl")

# Distributed order differential equations
include("DistributedOrder/utils.jl")
include("DistributedOrder/matrix.jl")

# Fractional Differences equations
include("FractionalDifferences/PECE.jl")

# Fractional integral equations
include("FIE/Qmat.jl")
include("FIE/mycoeffs.jl")
include("FIE/main.jl")

# Mittag Leffler function
include("mlfun.jl")


export FractionalDiffEqAlgorithm

# Export general types
export solve, FDEProblem, FPDEProblem, FDDEProblem

# Detailed problem types
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem, FractionalDifferenceProblem, FIEProblem, FDDESystem

# FODE solvers
export PIEx, PIIm, PITrap
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, GL
export IMEX_I, IMEX_I_A, IMEX_II

# FPDE solvers
export FPDEMatrixDiscrete, FiniteDiffEx, FiniteDiffIm, ADV_DIF, GLDiff

# System of FDE solvers
export NonLinearAlg, GLWithMemory, LorenzADM

# FDDE solvers
export DelayPECE, DelayPI, DelayABM, DelayABMYuan, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
export PECEDifference

# Fractional Integral Equations solvers
export SpectralUltraspherical



# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special models
export bagleytorvik, diffusion, FractionalLorenz

# Auxiliary functions
export mittleff, mittlefferr, mittleffderiv

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

export myeval

end