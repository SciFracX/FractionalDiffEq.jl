module FractionalDiffEq

using LinearAlgebra
using SpecialFunctions
using SparseArrays
using ApproxFun: Fun, Ultraspherical, coefficients
using InvertedIndices
using SpecialMatrices
using FFTW: fft, ifft
using UnPack: @unpack
using LoopVectorization: @turbo
using HypergeometricFunctions
using ToeplitzMatrices
using RecipesBase
using ForwardDiff

include("PECE.jl")
include("matrix.jl")
include("PI.jl")
include("SystemPI/PIEX.jl")
include("MTPI/MTPIEX.jl")
include("MTPI/MTPIPECE.jl")
include("MTPI/MTPIIMTrap.jl")
include("MTPI/MTPIIMRect.jl")
include("SystemPI/PIPECE.jl")
include("FLMM/FLMMBDF.jl")
include("FLMM/FLMMNewtonGregory.jl")
include("FLMM/FLMMTrap.jl")
include("ClosedForm/hankelmatrix.jl")
include("ClosedForm/ClosedForm.jl")
include("ClosedForm/highprecision.jl")

include("GL/GLWithMemory.jl")
include("GL/GL.jl")
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
include("FractionalDifferences/GL.jl")

# Fractional integral equations
include("FIE/Qmat.jl")
include("FIE/mycoeffs.jl")
include("FIE/main.jl")

# Mittag Leffler function
include("mlfun.jl")

include("utils.jl")
include("auxiliary.jl")


export FractionalDiffEqAlgorithm

# General types
export solve, FDEProblem

# FPDE problems
export FPDEProblem

# FDDE problems
export FDDEProblem, FDDESystem, FDDEMatrixProblem

# FODE problems
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem

# Fractional Difference probelms
export FractionalDifferenceProblem, FractionalDifferenceSystem

# FIE problems
export FIEProblem


###################################################


export AbstractFDESolution
export FODESolution, FIESolution, FDifferenceSolution, DODESolution

# FODE solvers
export PIEX, PIPECE, PIIMRect, PIIMTrap
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, ClosedFormHighPercision, GL
export ChebSpectral

# FPDE solvers
export FPDEMatrixDiscrete, FiniteDiffEx, FiniteDiffIm, ADV_DIF, GLDiff

# System of FODE solvers
export NonLinearAlg, GLWithMemory, FLMMBDF, FLMMNewtonGregory, FLMMTrap, PIEX

# FDDE solvers
export DelayPECE, DelayPI, DelayABM, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
export PECEDifference

# Fractional Integral Equations solvers
export SpectralUltraspherical

###################################################

# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special models
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittleffderiv

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

export myeval, ourfft, ourifft, rowfft, FastConv

end