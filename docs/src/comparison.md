# Comparison with other toolboxes

There are also some toolboxes can be used for fractional differential equations, here, we make a survey of all these toolboxes and have a more intuitive perspective.

## Matlab

### FOTF toolbox

FOTF toolbox is a Matlab toolbox developed by Prof Dingyu Xue, with fractional calculus, fractional differential equations and fractional order systems integrated together. To be honest, FOTF has a great impact on FractionalDiffEq.jl, all of the algorithm in FOTF toolbox is all supported in FractionalDiffEq.jl

#### corresponding APIs:

| FOTF toolbox | FractionalDiffEq.jl |
| -- | -- |
| ```fode_sol9``` | ```ClosedForm``` |
| | |

### Matrix approach to discretization of ODEs and PDEs of arbitrary real order

The corresponding matrix discretization paper and toolbox is developed by Prof Igor Podlubny. The corresponding **matrix discretization method** in FractionalDiffEq.jl is ```FODEMatrixDiscrete``` and ```FPDEMatrixDiscrete```.

## Python

Don't see any packages for fractional differential equations.

## R

Also don't see any packages for fractional differential equations.