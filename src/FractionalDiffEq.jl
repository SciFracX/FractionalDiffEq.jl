module FractionalDiffEq

using LinearAlgebra
using SciMLBase, DiffEqBase
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
using TruncatedStacktraces

include("types/problems.jl")
include("types/algorithms.jl")
include("types/solutions.jl")

include("types/problem_utils.jl")

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
include("fodesystem/atangana_seda.jl")

# System of fractal-fractional ordinary differential equations
include("ffode/atangana_seda.jl")

# Fractional delay differential equations
include("delay/DelayPECE.jl")
include("delay/product_integral.jl")
include("delay/matrix_form.jl")
include("delay/DelayABM.jl")

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
export FODEProblem, MultiTermsFODEProblem, DODEProblem, FFPODEProblem, FFEODEProblem, FFMODEProblem

# Fractional Discrete probelms
export FractionalDiscreteProblem, FractionalDiscreteSystem


###################################################


export AbstractFDESolution
export FODESolution, FDifferenceSolution, DODESolution, FFMODESolution
export FODESystemSolution, FDDESystemSolution, FFMODESystem

# FODE solvers
export PIPECE, PIRect, PITrap
export PECE, FODEMatrixDiscrete, ClosedForm, ClosedFormHankelM, GL
export AtanganaSedaAB
#export Euler

# System of FODE solvers
export NonLinearAlg, BDF, NewtonGregory, Trapezoidal, PIEX, NewtonPolynomial
export AtanganaSedaCF
export AtanganaSeda

# FDDE solvers
export DelayPECE, DelayABM, MatrixForm

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
export FOLyapunov, FOLE

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

end