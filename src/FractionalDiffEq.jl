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

# Single-term fractional ordinary differential equations
include("singletermfode/PECE.jl")
include("singletermfode/PI.jl")
include("singletermfode/GL.jl")
include("singletermfode/ChebSpectral.jl")
include("singletermfode/AtanganaSeda.jl")

# Multi-terms fractional ordinary differential equations
include("multitermsfode/matrix.jl")
include("multitermsfode/hankelmatrix.jl")
include("multitermsfode/ClosedForm.jl")
include("multitermsfode/highprecision.jl")
include("multitermsfode/MTPIEX.jl")
include("multitermsfode/MTPIPECE.jl")
include("multitermsfode/MTPIIMTrap.jl")
include("multitermsfode/MTPIIMRect.jl")

# System of fractional ordinary differential equations
include("fodesystem/PIPECE.jl")
include("fodesystem/FLMMBDF.jl")
include("fodesystem/FLMMNewtonGregory.jl")
include("fodesystem/FLMMTrap.jl")
include("fodesystem/PIEX.jl")
include("fodesystem/GLWithMemory.jl")
include("fodesystem/NonLinear.jl")
include("fodesystem/NewtonPolynomial.jl")
include("fodesystem/AS.jl")

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
include("dode/utils.jl")
include("dode/matrix.jl")

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
export AtanganaSeda, AtanganaSedaAB

# FPDE solvers
export FPDEMatrixDiscrete, FiniteDiffEx, FiniteDiffIm, ADV_DIF, GLDiff

# System of FODE solvers
export NonLinearAlg, GLWithMemory, FLMMBDF, FLMMNewtonGregory, FLMMTrap, PIEX, NewtonPolynomial

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