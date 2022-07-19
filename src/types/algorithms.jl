"""
Base type for all of the FractionalDiffEq algorithms
"""
abstract type AbstractFDEAlgorithm end

"""
Base type for distributed order differential equations algorithms.
"""
abstract type DODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional delay differential equations algorithms.
"""
abstract type FDDEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order integral equations algorithms.
"""
abstract type FIEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order partial differential equations algorithms.
"""
abstract type FPDEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for system of fractional order ordinary differential equations algorithms.
"""
abstract type FODESystemAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for multi-terms fractional ordinary differential equations algorithms.
"""
abstract type MultiTermsFODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for single-term fractional ordinary differential equations algorithms.
"""
abstract type SingleTermFODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order difference equations algorithms.
"""
abstract type FDifferenceAlgorithm <: AbstractFDEAlgorithm end