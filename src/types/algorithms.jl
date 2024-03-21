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
abstract type FODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for multi-terms fractional ordinary differential equations algorithms.
"""
abstract type MultiTermsFODEAlgorithm <: AbstractFDEAlgorithm end

"""
Base type for fractional order difference equations algorithms.
"""
abstract type FractionalDiscreteAlgorithm <: AbstractFDEAlgorithm end




###################### FODE ######################

"""
  BDF

[Backward Differentiation Formula](https://en.wikipedia.org/wiki/Backward_differentiation_formula) generated weights fractional linear multi-steps method.

## References

@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
"""
struct BDF <: FODEAlgorithm end

"""
Grunwald Letnikov difference method.

## Reference

@INPROCEEDINGS{8742063,  
author={Clemente-López, D. and Muñoz-Pacheco, J. M. and Félix-Beltrán, O. G. and Volos, C.},  
booktitle={2019 8th International Conference on Modern Circuits and Systems Technologies (MOCAST)},   
title={Efficient Computation of the Grünwald-Letnikov Method for ARM-Based Implementations of Fractional-Order Chaotic Systems},   
year={2019},   
doi={10.1109/MOCAST.2019.8742063}}

Python version by https://github.com/DClementeL/Grunwald_Letnikov
"""
# Grunwald Letnikov discretization method dispatch for FODEProblem
# struct GLWithMemory <: FractionalDiffEqAlgorithm end
struct GL <: FODEAlgorithm end


"""
    solve(prob::FODEProblem, NewtonGregory(); abstol=1e-3, maxiters=1e3)

[Newton Gregory](https://www.geeksforgeeks.org/newton-forward-backward-interpolation/) generated weights fractional linear multiple steps method.

## References

@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
"""
struct NewtonGregory <: FODEAlgorithm end


"""
    solve(prob::FODEProblem, dt, NewtonPolynomial())

Solve FODE system using Newton Polynomials methods.

!!! tip
    Used for the Caputo Fabrizio fractional differential operators.

## References

https://doi.org/10.1016/c2020-0-02711-8
"""
struct NewtonPolynomial <: FODEAlgorithm end


"""
# Usage

    solve(prob::FODEProblem, dt, NonLinearAlg())

Nonlinear algorithm for nonlinear fractional differential equations.

## References

Dingyu Xue, Northeastern University, China ISBN:9787030543981
"""
struct NonLinearAlg <: FODEAlgorithm end


"""
    solve(prob::FODEProblem, Trapezoidal(); abstol=1e-3, maxiters=1e3)

Use [Trapezoidal](https://en.wikipedia.org/wiki/Trapezoidal_rule_(differential_equations)) with generating function ``f(x)=\\frac{1+x}{2(1-x)^\\alpha}`` generated weights fractional linear multiple steps method to solve system of FODE.

## References

@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
"""
struct Trapezoid <: FODEAlgorithm end


"""
    Predictor-Correct scheme.

## References

@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
"""
struct PECE <: FODEAlgorithm end

"""
# Usage

    solve(FDDE::FDDEProblem, dt, DelayPECE())

Using the delayed predictor-corrector method to solve the delayed fractional differential equation problem in the Caputo sense.

Capable of solving both single term FDDE and multiple FDDE, support time varying lags of course😋.

## References

@article{Wang2013ANM,
  title={A Numerical Method for Delayed Fractional-Order Differential Equations},
  author={Zhen Wang},
  journal={J. Appl. Math.},
  year={2013},
  volume={2013},
  pages={256071:1-256071:7}
}

@inproceedings{Nagy2014NUMERICALSF,
  title={NUMERICAL SIMULATIONS FOR VARIABLE-ORDER FRACTIONAL NONLINEAR DELAY DIFFERENTIAL EQUATIONS},
  author={Abdelhameed M. Nagy and Taghreed Abdul Rahman Assiri},
  year={2014}
}

@inproceedings{Abdelmalek2019APM,
  title={A Predictor-Corrector Method for Fractional Delay-Differential System with Multiple Lags},
  author={Salem Abdelmalek and Redouane Douaifia},
  year={2019}
}
"""
struct DelayPECE <: FDDEAlgorithm end

"""
    Euler

The classical Euler method extended for fractional order differential equations.
"""
struct Euler <: FODEAlgorithm end

"""
  PIEX

Explicit product integral method for initial value problems of fractional order differential equations.
"""
struct PIEX <: FODEAlgorithm end

"""
  DelayPIEX

Explicit product integral method for initial value problems of fractional order differential equations.
"""
struct DelayPIEX <: FDDEAlgorithm end


"""
  MTPIEX

Explicit product integral method for initial value problems of fractional order differential equations.
"""
struct MTPIEX <: MultiTermsFODEAlgorithm end



"""
  AtanganaSedaAB

Solve Atangana-Baleanu fractional order differential equations using Newton Polynomials.
"""
struct AtanganaSedaAB <: FODEAlgorithm end


"""
  MatrixDiscrete

[Triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional ordinary differential equations to simple algebra system and solve the system.

## References

@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
"""
struct MatrixDiscrete <: MultiTermsFODEAlgorithm end