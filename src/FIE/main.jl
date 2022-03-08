"""
# Usage

    solve(f, e, n, SpectralUltraspherical())

### References

```tex
@article{Hale2018AFA,
  title={A Fast and Spectrally Convergent Algorithm for Rational-Order Fractional Integral and Differential Equations},
  author={Nicholas Hale and Sheehan Olver},
  journal={SIAM J. Sci. Comput.},
  year={2018},
  volume={40},
  pages={A2456-A2491}
}
```
"""
struct SpectralUltraspherical <: FractionalDiffEqAlgorithm end

function solve(f, e, n, ::SpectralUltraspherical)
    # Initialize operators
    Z = spzeros(n, n)
    I1 = sparse(Matrix(I, n, n))
    II = [I1 Z; Z I1]
    Q05(lam) = Qmat(n, 0.5, lam)
    QQ05 = [Z Q05(1); Q05(0.5) Z]
    QQ = m->QQ05^(2*m)

    A = II+QQ05
    rhs = [mycoeffs(e, n, 0.5); mycoeffs(f, n, 1)]

    idx = [collect(1:n)'; collect(n+1:2*n)']
    idx=idx[:]
    A = A[idx, idx]

    # Initialize solution
    u=zeros(size(A, 1))
    u[idx] = A\rhs
    return u
end