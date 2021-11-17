import FractionalDiffEq.solve, FractionalDiffEq.FractionalDiffEqAlgorithm

using LinearAlgebra, InvertedIndices

"""
Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional differential equations to simple algebra system and solve the system.

```tex
@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
```
"""
struct MatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
@inproceedings{Podlubny1998FractionalDE,
  title={Fractional differential equations},
  author={Igor Podlubny},
  year={1998}
}
"""

"""
@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
"""


"""
    solve(equation, right, h, T, MatrixDiscrete())

Using the **Matrix Discretization algorithm** proposed by [Prof Igor Podlubny](http://people.tuke.sk/igor.podlubny/index.html) to approximate the numerical solution.
"""
function solve(equation, right, h, T, ::MatrixDiscrete)
    N=Int64(T/h+1)
    equation = eliminator(N, rows)*equation*eliminator(N, rows)'
    right = right*eliminator(N, rows)*ones(N)
    result = equation\right
    return vcat(zeros(n), result)
end

"""
    eliminator(n, row)

Compute the eliminator matrix Sₖ by omiting n-th row
"""
function eliminator(n, row)
    temp = zeros(n, n)+I
    return temp[Not(row), :]
end

"""
Generating elements in Matrix.
"""
function omega(n, α)
    omega = zeros(n+1)

    omega[1]=1
    for i in range(1, n, step=1)
        omega[i+1]=(1-(α+1)/i)*omega[i]
    end
    
    return omega

end

"""
    D(N, α, h)

Using this function to construct left hand equations.

### Example

Suppose we have a equation:

```math
y''(t)+D^{frac{3}{2}}_t y(t)=1
```

```julia-repl
equation=D(100, 2, 0.01)+D(100, 3/2, 0.01)
```

We construct the left side equation!

After we construct the right side, alll we need to do is:

```julia-repl
solve(equation, right, h, T)
```

And then we can get the numerical approximation.

!!! info
    Here ```N``` is the size of discrete matrix.
"""
function D(N, α, h)
    result=zeros(N, N)
    temp=omega(N, α)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    return h^(-α)*result
end

function F(N, α, h)
    result=zeros(N, N)
    temp = omega(N, α)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    result=reverse(reverse(result, dims=1), dims=2)

    return (-1)^(ceil(α))*h^(-α)*result
end

"""
    bagleytorvik(p1, p2, p3, T, h)

By specifying the parameters of Bagley Torvik Equation, we can use **bagleytorvik** to directly obtain the numerical approximation.

!!! info "p2 ≠ 0"
    Please note that the parameter of fractional derivative item must not be 0
"""
function bagleytorvik(p1, p2, p3, T, h)
    N=Int64(T/h+1)
    equation = p1*D(N, 2, h)+p2*D(N, 1.5, h)+p3*(zeros(N, N)+I)
    equation = eliminator(N, [1,2])*equation*eliminator(N, [1,2])'
    right = eliminator(N, [1,2])*ones(N)
    result = equation\right

    return vcat(zeros(2), result)
end