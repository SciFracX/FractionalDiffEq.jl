# FractionalDiffEq.jl APIs

## Problem definition:

```@docs
FractionalDiffEq.FDEProblem
```

## General solve API

```@docs
FractionalDiffEq.solve
```

The general solving API ```solve``` can accept various kinds of inputs, including [FDEProblem](@ref FDEProblem), [FODEProblem](@ref FODEProblem) and [FPDEProblem](@ref FPDEProblem)

## Current algorithms:

### Base type

```@docs
FractionalDiffEq.FractionalDiffEqAlgorithm
```

## Problem types

```@docs
FractionalDiffEq.FODEProblem
FractionalDiffEq.FPDEProblem
```

### Detailed types

```@docs
FractionalDiffEq.PECE
```

```@docs
FractionalDiffEq.FODEMatrixDiscrete
FractionalDiffEq.FPDEMatrixDiscrete
```

## Some models

```@docs
FractionalDiffEq.bagleytorvik
```

```@docs
FractionalDiffEq.diffusion
```