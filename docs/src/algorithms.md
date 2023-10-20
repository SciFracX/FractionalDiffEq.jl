# Current algorithms in FractionalDiffEq.jl

```@contents
Pages = ["algorithms.md"]
```

## Fractional Ordinary Differential Equations

### Single Term FODE

```@docs
FractionalDiffEq.PECE
FractionalDiffEq.Euler
FractionalDiffEq.PIEX
FractionalDiffEq.AtanganaSeda
```

### Multi-Term FODE

```@docs
FractionalDiffEq.FODEMatrixDiscrete
FractionalDiffEq.ClosedForm
FractionalDiffEq.ClosedFormHankelM
FractionalDiffEq.PIPECE
FractionalDiffEq.PIRect
FractionalDiffEq.PITrap
```

### System of FODE

```@docs
FractionalDiffEq.NonLinearAlg
FractionalDiffEq.PECE
FractionalDiffEq.GL
FractionalDiffEq.FLMMBDF
FractionalDiffEq.FLMMNewtonGregory
FractionalDiffEq.FLMMTrap
FractionalDiffEq.PIEX
FractionalDiffEq.NewtonPolynomial
FractionalDiffEq.AtanganaSedaAB
```

## Fractional Delay Differential Equatinos

```@docs
FractionalDiffEq.DelayPECE
FractionalDiffEq.DelayPI
FractionalDiffEq.MatrixForm
FractionalDiffEq.DelayABM
```

## Distributed Order Differential Equations

```@docs
FractionalDiffEq.DOMatrixDiscrete
```

## Fractional Differences Equations

```@docs
FractionalDiffEq.PECEDifference
FractionalDiffEq.GL
```