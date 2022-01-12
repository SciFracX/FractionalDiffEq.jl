```@meta
CurrentModule = FractionalDiffEq
DocTestSetup  = quote
    using FractionalDiffEq
end
DocTestFilters = [r"Stacktrace:[\s\S]+"]
```

# FractionalDiffEq.jl APIs

## Problem definition:

```@docs
FractionalDiffEq.FDEProblem
```

## General solve API

```@docs
FractionalDiffEq.solve
```

The general solving API [```solve```](#FractionalDiffEq.solve) can accept various kinds of inputs, including [FDEProblem](#FractionalDiffEq.FDEProblem), [FODEProblem](#FractionalDiffEq.FODEProblem) and [FPDEProblem](#FractionalDiffEq.FPDEProblem)

## Current algorithms:

### Base type

```@docs
FractionalDiffEq.FractionalDiffEqAlgorithm
```

## Problem types

```@docs
FractionalDiffEq.SingleTermFODEProblem
FractionalDiffEq.MultiTermsFODEProblem
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