abstract type AbstractFDESolution end

struct FODESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FODESystemSolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FDDESystemSolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FDifferenceSolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FIESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct DODESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FFMODESolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end

struct FFMODESystemSolution <: AbstractFDESolution
    t::AbstractArray
    u::AbstractArray
end




struct FOLE
    t::AbstractArray
    LE::AbstractArray
end