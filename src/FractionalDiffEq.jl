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

include("problems.jl")
include("algorithms.jl")
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
export FDDEProblem, FDDEMatrixProblem

# FODE problems
export FODEProblem, MultiTermsFODEProblem, DODEProblem

# Fractional Discrete probelms
export FractionalDiscreteProblem, FractionalDiscreteSystem

###################################################

# FODE solvers
export PIRect, PITrap, PECE, PIEX
export MatrixDiscrete, GL

export MTPIRect, MTPITrap, MTPECE, MTPIEX

# System of FODE solvers
export NonLinearAlg, BDF, NewtonGregory, Trapezoid, NewtonPolynomial

# FDDE solvers
export DelayPIEX, DelayPECE, DelayABM, MatrixForm

# DODE solvers
export DOMatrixDiscrete

# export extension solvers
export FdeSolverPECE

###################################################

# Export some special models
export bagleytorvik, diffusion

# Auxiliary functions
export mittleff, mittleffderiv

# Lyapunov exponents
#export FOLyapunov, FOLE

# Distributed order auxiliary SpecialFunctions
export DOB, DOF, DORANORT, isFunction

#=
function chua!(du, x, p, t)
    a, b, c, m0, m1 = p
    du[1] = a*(x[2]-x[1]-(m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))))
    du[2] = x[1]-x[2]+x[3]
    du[3] = -b*x[2]-c*x[3]
end
α = [0.93, 0.99, 0.92];
x0 = [0.2, -0.1, 0.1];
tspan = (0, 100);
p = [10.725, 10.593, 0.268, -1.1726, -0.7872]
prob = FODEProblem(chua!, α, x0, tspan, p)
=#

end
