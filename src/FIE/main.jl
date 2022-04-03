"""
# Usage

    solve(prob::FIEProblem, n, SpectralUltraspherical())

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

function solve(FIE::FIEProblem, n, ::SpectralUltraspherical)
  @unpack parameters, orders, rightfun, tspan = FIE
  # Initialize operators
  Z = spzeros(n, n)
  I1 = sparse(Matrix(I, n, n))
  II = [I1 Z; Z I1]
  Q05(lam) = Qmat(n, 0.5, lam)
  QQ05 = [Z Q05(1); Q05(0.5) Z]
  QQ = m->QQ05^(2*m)

  judgeorders(j) = j == 1 ? (return II) : (return QQ(j))

  A = zeros(2*n, 2*n)

  for (i, j) in zip(parameters, orders)
    A += i*judgeorders(j)
  end

  f = 0
  rhs = [mycoeffs(rightfun, n, 0.5); mycoeffs(f, n, 1)]

  idx = [collect(1:n)'; collect(n+1:2*n)']
  idx = idx[:]
  A = A[idx, idx]

  # Initialize solution
  u=zeros(size(A, 1))
  u[idx] = A\rhs

  # Evaluate solution
  sol = myeval(u, tspan, 2)
  return sol
end

function myeval(u, x, q)
  if size(u, 2) > 1
      uu = zeros(length(x), size(u, 2))
      for k = 1:size(u, 2)
          uu[:, k] = myeval(u[:, k], x, q)
      end
      return uu
  end
  
  N::Int64 = length(u)/q
  
  if q == 2
      uu = clenshawP(x, u[1:N]) .+ sqrt.(1 .+x).*clenshawU(x, u[N+1:2*N])
  else
      uu = 0*x
      for k = 0:q-1
          uu = uu .+ (1+x).^(k/q).*clenshawJ_special(x, u[k*N+collect(1:N)], k/q)
      end
  end
  return uu
end
  
function clenshawP(x, c)
  bk1 = 0*x
  bk2 = copy(bk1)
  N = size(c, 1)-1
  for k = N:-1:1
      bk = c[k+1] .+ (2*k+1)/(k+1)*x.*bk1 .- (k+1)/(k+2)*bk2
      bk2 = bk1
      bk1 = bk
  end
  y = c[1] .+ x.*bk1 - 0.5*bk2
  return y
end
  
function clenshawU(x, c)
  bk1 = 0*x
  bk2 = copy(bk1)
  N = size(c, 1)-1
  for k = N:-1:1
      bk = c[k+1] .+ 2*x.*bk1 - bk2
      bk2 = bk1
      bk1 = bk
  end
  y = c[1] .+ 2*x.*bk1 - bk2
  return y
end
  
function clenshawJ_special(x, c, b)
  bk1 = 0*x
  bk2 = copy(bk1)
  N = size(c, 1)-1
  for k = N:-1:1
      Ak = (2*k+3)
      Bk = (1-2*b)/(2*k+1)
      Ck = (k+2-b)*(k+1+b)*(2*k+5)/((k+3)*(2*k+3))
      bk = c[k+1] + ((Ak*x+Bk).*bk1 - Ck*bk2)/(k+2)
      bk2 = bk1
      bk1 = bk
  end
  Ck = (2-b)*(1+b)*5/18
  y = c[1] + 0.5*(3*x+1-2*b).*bk1 - Ck*bk2
  return y
end