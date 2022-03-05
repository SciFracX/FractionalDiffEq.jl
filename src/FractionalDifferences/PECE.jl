"""
# Usage

    solve(FDProb::FractionalDifferenceProblem, T, h, PECEDifference())

Use the PECE algorithm to solve fractional difference equations

### References

```tex
@article{陈福来2011分数阶微分差分方程的,
  title={分数阶微分差分方程的 Matlab 求解},
  author={陈福来 and 王华 and 李势丰},
  journal={湘南学院学报},
  volume={32},
  number={5},
  pages={1--4},
  year={2011}
}
```
"""
struct PECEDifference <: FractionalDiffEqAlgorithm end


function solve(FDProb::FractionalDifferenceProblem, T, h, ::PECEDifference)
    @unpack fun, α, x0 = FDProb
    N = round(Int, T/h)
    b=zeros(N); a=zeros(N); x=zeros(N); t=zeros(N)
    t1=zeros(N+1); y=zeros(N+1)
    @fastmath @inbounds @simd for i=1:N
        k = 1
        @fastmath @inbounds @simd for j=1:N-i
            k = k*j
        end
        b[i] = gamma(N-i+α)/k
        a[N+1-i] = b[i]
    end
    x[1] = x0+a[1]*fun(x0)/gamma(α)
    t[1] = 1

    @fastmath @inbounds @simd for i=2:N
        temp = 0
        @fastmath @inbounds @simd for j=1:i-1
            temp = a[i+1-j]*fun(x[j])
        end
        temp = temp/gamma(α)
        x[i] = x0 + temp+a[1]*(x0+temp)/gamma(α)
        t[i] = i
    end

    y[1] = x0
    t1[1] = 0
    @fastmath @inbounds @simd for i = 1:N
        y[i+1] = x[i]
        t1[i+1] = t[i]
    end
    return t1, y
end