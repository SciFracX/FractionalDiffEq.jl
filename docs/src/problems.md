# Problem Types

```@contents
Pages = ["problems.md"]
```

The general fractional differential equations problem type:

```@docs
FractionalDiffEq.FDEProblem
```

## FODE Problem

Fractional ordinary problems type, we can classified them as Single-Term and Multi-Term problems:

```math
D^{\alpha}y(t)=f(t, y)
```

```math
\sum_{s=0}^pc_sD^{\beta_s}y=f
```

```@docs
FractionalDiffEq.SingleTermFODEProblem
FractionalDiffEq.MultiTermsFODEProblem
FractionalDiffEq.FODESystem
```

## FPDE Problem

```@docs
FractionalDiffEq.FPDEProblem
```

## FDDEProblem

```@docs
FractionalDiffEq.FDDEProblem
FractionalDiffEq.FDDEMatrixProblem
FractionalDiffEq.FDDESystem
```

## DODEProblem

```@docs
FractionalDiffEq.DODEProblem
```

## Fractional Difference Problem

```@docs
FractionalDiffEq.FractionalDifferenceProblem
FractionalDiffEq.FractionalDifferenceSystem
```

## FIEProblem

```@docs
FractionalDiffEq.FIEProblem
```