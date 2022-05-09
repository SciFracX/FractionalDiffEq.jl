module FractionalDiffEq

using LinearAlgebra, SpecialFunctions, SparseArrays, ApproxFun, InvertedIndices, QuadGK
using SpecialMatrices, Polynomials, FFTW, UnPack, LoopVectorization, HypergeometricFunctions
using ToeplitzMatrices, RecipesBase, ForwardDiff

include("PECE.jl")
include("matrix.jl")
include("PI.jl")
include("PIEX.jl")
include("MTPI/MTPIEX.jl")
include("MTPI/MTPIPECE.jl")
include("MTPI/MTPIIMTrap.jl")
include("MTPI/MTPIIMRect.jl")
include("PIPECE.jl")
include("FLMM/FLMMBDF.jl")
include("FLMM/FLMMNewtonGregory.jl")
include("FLMM/FLMMTrap.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

include("GL/GLWithMemory.jl")
include("GL/GL.jl")
#include("Direct.jl")
include("ChebSpectral.jl")

include("NonLinear/NonLinear.jl")

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
include("FDDE/DelayABMSystem.jl")

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

include("utils.jl")
include("auxiliary.jl")


export FractionalDiffEqAlgorithm

# Export general types
export solve, FDEProblem, FPDEProblem, FDDEProblem

# Detailed problem types
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem, FractionalDifferenceProblem, FIEProblem, FDDESystem, FDDEMatrixProblem

export AbstractFDESolution, FODESolution, FIESolution, FDifferenceSolution, DODESolution

# FODE solvers
export PIEX, PIIM, PITrap, PIPECE, PIIMRect, PIIMTrap
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, GL
export ChebSpectral

# FPDE solvers
export FPDEMatrixDiscrete, FiniteDiffEx, FiniteDiffIm, ADV_DIF, GLDiff

# System of FDE solvers
export NonLinearAlg, GLWithMemory, FLMMBDF, FLMMNewtonGregory, FLMMTrap, PIEX

# FDDE solvers
export DelayPECE, DelayPI, DelayABM, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
export PECEDifference

# Fractional Integral Equations solvers
export SpectralUltraspherical



# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special models
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittlefferr, mittleffderiv

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

export myeval

end