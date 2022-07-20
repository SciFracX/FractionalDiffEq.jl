"""
### References

```tex
@article{Wei2022HowTE,
  title={How to empower Gr{\"u}nwald–Letnikov fractional difference equations with available initial condition?},
  author={Yiheng Wei and Jinde Cao and Chuang Li and Yang Quan Chen},
  journal={Nonlinear Analysis: Modelling and Control},
  year={2022}
}
```
"""


function solve(prob::FractionalDifferenceSystem, N, ::GL)
    @unpack fun, α, u0, p = prob
    result = zeros(Float64, length(u0), N)# Initialization
    result[:, 1] = u0

    for j=2:N
        du = zeros(Float64, length(u0))
        fun(du, result[:, j-1], p, 0)
        result[:, j] = du
        # I think this way would be faster, but if I pass the result[:, j] to fun,
        # that would create a copy, which is not what I want.
        # fun(result[:, j], result[:, j-1], nothing, 0)
        for i=1:j-2
            result[:, j] = result[:, j] - generalized_binomial(i, α).*result[:, j-i]
        end
    end
    return result
end


function generalized_binomial(i, alpha)
    result = ones(Float64, i+1)
    for j in 1:i
        result[j+1] = result[j]*(j-1-alpha)/j
    end
    return result[end]
end