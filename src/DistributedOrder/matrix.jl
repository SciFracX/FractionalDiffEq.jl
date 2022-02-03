import FractionalDiffEq: FractionalDiffEqAlgorithm, solve, eliminator, DOB

using LinearAlgebra

"""

### References

```tex
@inproceedings{Jiao2012DistributedOrderDS,
  title={Distributed-Order Dynamic Systems - Stability, Simulation, Applications and Perspectives},
  author={Zhuang Jiao and Yang Quan Chen and Igor Podlubny},
  booktitle={Springer Briefs in Electrical and Computer Engineering},
  year={2012}
}
```
"""
struct DOMatrixDiscrete <: FractionalDiffEqAlgorithm end

function solve(ω, t, h, B, ::DOMatrixDiscrete)
    N = length(t)+1
    M = zeros(N, N)

    M = DOB(ω, [0, 1], 0.01, N-1, h)+B*(zeros(N-1, N-1)+I)
    F = -ones(N-1)./10

    M = eliminator(N-1, 1)*M*eliminator(N-1, 1)'
    F = eliminator(N-1, 1)*F

    Y = M\F

    Y0 = [0; Y]

    return Y0.+1
end
