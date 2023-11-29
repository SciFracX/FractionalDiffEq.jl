# Problem Types

```@contents
Pages = ["problems.md"]
```

The general fractional differential equations problem type.

## FODE Problem

Fractional ordinary problems type, we can classify them as Single-Term and Multi-Term problems:

```math
D^{\alpha}y(t)=f(t, y)
```

```math
\sum_{s=0}^pc_sD^{\beta_s}y=f
```

SingleTermFODEProblem
MultiTermsFODEProblem
FODESystem

### FFODE Problem

FFPODEProblem
FFEODEProblem
FFMODEProblem

## FDDE Problem

FDDEProblem
FDDEMatrixProblem
FDDESystem

## DODE Problem

FractionalDiffEq.DODEProblem

## Fractional Difference Problem

FractionalDifferenceProblem
FractionalDifferenceSystem

## FIE Problem

FIEProblem
