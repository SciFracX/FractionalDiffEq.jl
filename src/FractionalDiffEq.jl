module FractionalDiffEq

using LinearAlgebra, Reexport, SciMLBase, SpecialFunctions, SparseArrays, ToeplitzMatrices, RecursiveArrayTools,
      FFTW, ForwardDiff, Polynomials, HypergeometricFunctions, DiffEqBase, ConcreteStructs, FastClosures

import SciMLBase: __solve
import DiffEqBase: solve
import InvertedIndices: Not
import SpecialMatrices: Vandermonde
import FFTW: fft, ifft
import Polynomials: Polynomial

@reexport using DiffEqBase, SciMLBase

include("types/problems.jl")
include("types/algorithms.jl")
include("types/solutions.jl")

# Multi-terms fractional ordinary differential equations
include("multitermsfode/matrix.jl")
include("multitermsfode/hankelmatrix.jl")
include("multitermsfode/closed_form.jl")
include("multitermsfode/explicit_pi.jl")
include("multitermsfode/implicit_pi_pece.jl")
include("multitermsfode/implicit_pi_trapezoid.jl")
include("multitermsfode/implicit_pi_rectangle.jl")

# System of fractional ordinary differential equations
include("fode/pi_pece.jl")
include("fode/bdf.jl")
include("fode/newton_gregory.jl")
include("fode/trapezoid.jl")
include("fode/explicit_pi.jl")
include("fode/implicit_pi_rectangle.jl")
include("fode/implicit_pi_trapzoid.jl")
include("fode/grunwald_letnikov.jl")
include("fode/nonlinearalg.jl")

# System of fractal-fractional ordinary differential equations
include("ffode/atangana_seda.jl")

# Fractional delay differential equations
include("delay/pece.jl")
include("delay/product_integral.jl")
include("delay/matrix_form.jl")
include("delay/abm.jl")

# Distributed order differential equations
include("dode/utils.jl")
include("dode/matrix.jl")

# Fractional Differences equations
include("discrete/PECE.jl")
include("discrete/GL.jl")

# Mittag Leffler function
include("mlfun.jl")

# Lyapunov exponents
#include("lyapunov.jl")

include("utils.jl")
include("auxiliary.jl")

function __solve(
        prob::MultiTermsFODEProblem, alg::MultiTermsFODEAlgorithm, args...; kwargs...)
    cache = init(prob, alg, args...; kwargs...)
    return solve!(cache)
end

function __solve(prob::FODEProblem, alg::FODEAlgorithm, args...; kwargs...)
    cache = init(prob, alg, args...; kwargs...)
    return solve!(cache)
end

function __solve(prob::FDDEProblem, alg::FDDEAlgorithm, args...; kwargs...)
    cache = init(prob, alg, args...; kwargs...)
    return solve!(cache)
end

# General types
export FDEProblem

# FDDE problems
export FDDEProblem, FDDESystem, FDDEMatrixProblem

# FODE problems
export FODEProblem, MultiTermsFODEProblem, DODEProblem, FFPODEProblem, FFEODEProblem,
       FFMODEProblem

# Fractional Discrete probelms
export FractionalDiscreteProblem, FractionalDiscreteSystem

###################################################

export AbstractFDESolution
export FODESolution, FDifferenceSolution, DODESolution, FFMODESolution
export FODESystemSolution, FDDESystemSolution, FFMODESystem

# FODE solvers
export PIRect, PITrap, PECE, PIEX
export MatrixDiscrete, GL
export AtanganaSedaAB

export MTPIRect, MTPITrap, MTPECE, MTPIEX

# System of FODE solvers
export NonLinearAlg, BDF, NewtonGregory, Trapezoid, NewtonPolynomial

# FDDE solvers
export DelayPIEX, DelayPECE, DelayABM, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# Fractional Differences Equations solvers
# export PECE

###################################################

# Export some special models
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittleffderiv

# Lyapunov exponents
#export FOLyapunov, FOLE

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

end
