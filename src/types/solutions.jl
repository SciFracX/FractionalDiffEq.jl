abstract type AbstractFDESolution end

struct FODESolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FODESystemSolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FDDESystemSolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FDifferenceSolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FIESolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct DODESolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FFMODESolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end

struct FFMODESystemSolution <: AbstractFDESolution
    t::AbstractMatrix
    u::AbstractMatrix
end