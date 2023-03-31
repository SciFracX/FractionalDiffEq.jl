module FractionalDiffEq

using LinearAlgebra
using SpecialFunctions
using SparseArrays
using InvertedIndices: Not
using SpecialMatrices: Vandermonde
using FFTW: fft, ifft
using UnPack: @unpack
using LoopVectorization: @turbo
using HypergeometricFunctions
using ToeplitzMatrices
using RecipesBase
using ForwardDiff
using Polynomials: Polynomial

include("types/problems.jl")
include("types/algorithms.jl")
include("types/solutions.jl")


# Single-term fractional ordinary differential equations
include("singletermfode/PECE.jl")
include("singletermfode/PI.jl")
include("singletermfode/GL.jl")
include("singletermfode/AtanganaSeda.jl")
include("singletermfode/Euler.jl")

# Multi-terms fractional ordinary differential equations
include("multitermsfode/matrix.jl")
include("multitermsfode/hankelmatrix.jl")
include("multitermsfode/closed_form.jl")
include("multitermsfode/explicit_pi.jl")
include("multitermsfode/implicit_pi_pece.jl")
include("multitermsfode/implicit_pi_trapezoid.jl")
include("multitermsfode/implicit_pi_rectangle.jl")

# System of fractional ordinary differential equations
include("fodesystem/PIPECE.jl")
include("fodesystem/bdf.jl")
include("fodesystem/newton_gregory.jl")
include("fodesystem/trapezoid.jl")
include("fodesystem/explicit_pi.jl")
include("fodesystem/Euler.jl")
include("fodesystem/GLWithMemory.jl")
include("fodesystem/NonLinear.jl")
include("fodesystem/newton_polynomials.jl")
include("fodesystem/AS.jl")
include("fodesystem/ASCF.jl")

# System of fractal-fractional ordinary differential equations
include("ffode/AS.jl")

# Fractional delay differential equations
include("delay/DelayPECE.jl")
include("delay/DelayPI.jl")
include("delay/Matrix.jl")
include("delay/DelayABM.jl")
include("delay/DelayABMSystem.jl")

# Distributed order differential equations
include("dode/utils.jl")
include("dode/matrix.jl")

# Fractional Differences equations
include("discrete/PECE.jl")
include("discrete/GL.jl")

# Mittag Leffler function
include("mlfun.jl")

# Lyapunov exponents
include("lyapunov.jl")

include("utils.jl")
include("auxiliary.jl")

# General types
export solve, FDEProblem

# FDDE problems
export FDDEProblem, FDDESystem, FDDEMatrixProblem

# FODE problems
export SingleTermFODEProblem, MultiTermsFODEProblem, FODESystem, DODEProblem, FFPODEProblem, FFEODEProblem, FFMODEProblem

# Fractional Discrete probelms
export FractionalDiscreteProblem, FractionalDiscreteSystem


###################################################


export AbstractFDESolution
export FODESolution, FDifferenceSolution, DODESolution, FFMODESolution
export FODESystemSolution, FDDESystemSolution, FFMODESystem

# FODE solvers
export PIEX, PIPECE, PIIMRect, PIIMTrap
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, GL
export AtanganaSeda, AtanganaSedaAB
export Euler

# System of FODE solvers
export NonLinearAlg, FLMMBDF, FLMMNewtonGregory, FLMMTrap, PIEX, NewtonPolynomial
export AtanganaSedaCF

# FDDE solvers
export DelayPECE, DelayPI, DelayABM, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
# export PECE

###################################################

# Export some api to construct the equation
export eliminator, RieszMatrix, omega, meshgrid

# Export some special models
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittleffderiv

# Lyapunov exponents
export FOLyapunov

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

export ourfft, ourifft, rowfft, FastConv

end