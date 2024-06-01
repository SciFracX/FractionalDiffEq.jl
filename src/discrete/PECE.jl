#=
# Usage

    solve(FDProb::FractionalDiscreteProblem, T, h, PECE())

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
=#

function solve(FDProb::FractionalDiscreteProblem, h, ::PECE)
    @unpack fun, α, u0, tspan = FDProb
    if tspan isa Tuple
        T = tspan[end]
    elseif tspan isa Real
        T = copy(tspan)
    end
    N = round(Int, T / h)

    # Initialize
    b = zeros(N)
    a = zeros(N)
    x = zeros(N)
    t = zeros(N)
    t1 = zeros(N + 1)
    y = zeros(N + 1)

    @fastmath @inbounds @simd for i in 1:N
        k = 1
        @fastmath @inbounds @simd for j in 1:(N - i)
            k = k * j
        end
        b[i] = gamma(N - i + α) / k
        a[N + 1 - i] = b[i]
    end
    x[1] = u0 + a[1] * fun(u0) / gamma(α)
    t[1] = 1

    @fastmath @inbounds @simd for i in 2:N
        temp = 0
        @fastmath @inbounds @simd for j in 1:(i - 1)
            temp = a[i + 1 - j] * fun(x[j])
        end
        temp = temp / gamma(α)
        x[i] = u0 + temp + a[1] * (u0 + temp) / gamma(α)
        t[i] = i
    end

    y[1] = u0
    t1[1] = 0
    @fastmath @inbounds @simd for i in 1:N
        y[i + 1] = x[i]
        t1[i + 1] = t[i]
    end
    return FDifferenceSolution(t1, y)
end
