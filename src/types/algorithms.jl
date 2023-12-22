"""
Base type for all of the FractionalDiffEq algorithms
"""
abstract type AbstractFDEAlgorithm <: SciMLBase.AbstractDEAlgorithm end

"""
Base type for distributed order differential equations algorithms.
"""
abstract type DODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional delay differential equations algorithms.
"""
abstract type FDDEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order partial differential equations algorithms.
"""
abstract type FPDEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for system of fractional order ordinary differential equations algorithms.
"""
abstract type FODESystemAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for multi-terms fractional ordinary differential equations algorithms.
"""
abstract type MultiTermsFODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for single-term fractional ordinary differential equations algorithms.
"""
abstract type SingleTermFODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order difference equations algorithms.
"""
abstract type FractionalDiscreteAlgorithm <: AbstractFDEAlgorithm end




###################### FODE ######################

"""
    solve(prob::FODEProblem, dt, BDF())

Use [Backward Differentiation Formula](https://en.wikipedia.org/wiki/Backward_differentiation_formula) generated weights fractional linear multi-steps method to solve system of FODE.

### References

```tex
@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
```
"""
struct BDF <: FODESystemAlgorithm end

"""
# Usage

    solve(prob::FODEProblem, dt, GL())

Use Grunwald Letnikov difference method to solve system of system of FODE.

### Reference

```tex
@INPROCEEDINGS{8742063,  
author={Clemente-López, D. and Muñoz-Pacheco, J. M. and Félix-Beltrán, O. G. and Volos, C.},  
booktitle={2019 8th International Conference on Modern Circuits and Systems Technologies (MOCAST)},   
title={Efficient Computation of the Grünwald-Letnikov Method for ARM-Based Implementations of Fractional-Order Chaotic Systems},   
year={2019},   
doi={10.1109/MOCAST.2019.8742063}}
```

Python version by https://github.com/DClementeL/Grunwald_Letnikov
"""
# Grunwald Letnikov discretization method dispatch for FODEProblem
# struct GLWithMemory <: FractionalDiffEqAlgorithm end
struct GL <: FODESystemAlgorithm end


"""
    solve(prob::FODEProblem, dt, FLMMNewtonGregory())

Use [Newton Gregory](https://www.geeksforgeeks.org/newton-forward-backward-interpolation/) generated weights fractional linear multiple steps method to solve system of FODE.

### References

```tex
@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
```
"""
struct NewtonGregory <: FODESystemAlgorithm end


"""
    solve(prob::FODEProblem, dt, NewtonPolynomial())

Solve FODE system using Newton Polynomials methods.

!!! tip
    Used for the Caputo Fabrizio fractional differential operators.

```tex
https://doi.org/10.1016/c2020-0-02711-8
```
"""
struct NewtonPolynomial <: FODESystemAlgorithm end


"""
# Usage

    solve(prob::FODEProblem, dt, NonLinearAlg())

Nonlinear algorithm for nonlinear fractional differential equations.

### References

Dingyu Xue, Northeastern University, China ISBN:9787030543981
"""
struct NonLinearAlg <: FODESystemAlgorithm end


"""
    solve(prob::FODEProblem, Trapezoidal())

Use [Trapezoidal](https://en.wikipedia.org/wiki/Trapezoidal_rule_(differential_equations)) with generating function ``f(x)=\\frac{1+x}{2(1-x)^\\alpha}`` generated weights fractional linear multiple steps method to solve system of FODE.

### References

```tex
@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
```
"""
struct Trapezoid <: FODESystemAlgorithm end


"""
Predictor-Correct scheme.

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
struct PECE <: FODESystemAlgorithm end


"""
The classical Euler method extended for fractional order differential equations.
"""
struct Euler <: FODESystemAlgorithm end



"""
    solve(prob::FODESystem, h, AtanganaSedaAB())

Solve Atangana-Baleanu fractional order differential equations using Newton Polynomials.
"""
struct AtanganaSedaAB <: FODESystemAlgorithm end


###################### FDDE ######################

"""
    solve(FDDE::FDDEProblem, dt, DelayABM())

Use the Adams-Bashforth-Moulton method to solve fractional delayed differential equations.

### References

```tex
@inproceedings{Bhalekar2011APS,
  title={A PREDICTOR-CORRECTOR SCHEME FOR SOLVING NONLINEAR DELAY DIFFERENTIAL EQUATIONS OF FRACTIONAL ORDER},
  author={Sachin Bhalekar and Varsha Daftardar-Gejji},
  year={2011}
}
```
"""
struct DelayABM <: FDDEAlgorithm end