# Comparison with other toolboxes

There are also some toolboxes that can be used for fractional differential equations, here we make a survey of all these toolboxes and have a more intuitive perspective.

## Matlab

### FOTF toolbox

FOTF toolbox is a Matlab toolbox developed by Prof Dingyu Xue, with fractional calculus, fractional differential equations and fractional order systems integrated together. To be honest, FOTF has a great impact on FractionalDiffEq.jl, all of the algorithms in FOTF toolbox is all supported in FractionalDiffEq.jl

#### corresponding APIs:

| FOTF toolbox | FractionalDiffEq.jl |
|-----|-----|
| ```fode_sol9``` | ```ClosedForm``` |
| ```fode_solm``` | ```ClosedFormHankelM``` |
| ```nlfode_vec``` | ```NonLinearAlg```|

### Matrix approach to discretization of ODEs and PDEs of arbitrary real order

The corresponding matrix discretization paper and toolbox is developed by Prof Igor Podlubny. The corresponding **matrix discretization method** in FractionalDiffEq.jl is ```FODEMatrixDiscrete``` and ```FPDEMatrixDiscrete```.

### FLMM2 Toolbox

FLMM2 toolbox is a toolbox developed by [Prof Roberto Garrappa](https://www.dm.uniba.it/members/garrappa) together with his [paper](http://dx.doi.org/10.1016/j.matcom.2013.09.012)

| Exsiting | FractionalDiffEq.jl |
|-----|-----|
| ```FLMM2(method=1)``` | ```FLMMTrap``` |
| ```FLMM2(method=2)``` | ```FLMMNewtonGregory``` |
| ```FLMM2(method=3)``` | ```FLMMBDF```|
| ```fde12``` | ```ABM``` |
| ```mt_fde_pi1_ex``` | ```PIEX``` |
| ```mt_fde_pi1_im``` | ```PIIMRect``` |
| ```mt_fde_pi2_im``` | ```PIIMTrap``` |
| ```mt_fde_pi12_pc``` | ```PIPECE``` |
## Python

[fodeint](https://github.com/mattja/fodeint): With explicit one-step Adams-Bashforth (Euler) method.

## R

Don't see any packages for fractional differential equations.