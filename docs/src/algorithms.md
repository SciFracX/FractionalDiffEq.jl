# Current algorithms in FractionalDiffEq.jl

```@contents
Pages = ["algorithms.md"]
```

## Fractional Ordinary Differential Equations

### Single Term FODE

```@docs
FractionalDiffEq.PECE
FractionalDiffEq.GL
FractionalDiffEq.PIEX
FractionalDiffEq.ChebSpectral
```

### Multi-Term FODE

```@docs
FractionalDiffEq.FODEMatrixDiscrete
FractionalDiffEq.ClosedForm
FractionalDiffEq.ClosedFormHankelM
FractionalDiffEq.PITrap
FractionalDiffEq.PIPECE
FractionalDiffEq.PIIMRect
FractionalDiffEq.PIIMTrap
```

### System of FODE

```@docs
FractionalDiffEq.NonLinearAlg
FractionalDiffEq.ABM
FractionalDiffEq.GL
FractionalDiffEq.FLMMBDF
FractionalDiffEq.FLMMNewtonGregory
FractionalDiffEq.FLMMTrap
FractionalDiffEq.PIEX
```

## Fractional Partial Differential Equations

```@docs
FractionalDiffEq.FPDEMatrixDiscrete
FractionalDiffEq.FiniteDiffEx
FractionalDiffEq.FiniteDiffIm
FractionalDiffEq.ADV_DIF
FractionalDiffEq.GLDiff
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

## Fractional Integral Equations

```@docs
FractionalDiffEq.SpectralUltraspherical
```