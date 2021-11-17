import FractionalDiffEq: solve, FractionalDiffEqAlgorithm

using LinearAlgebra, InvertedIndices

"""
@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
"""

"""
Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional ordinary differential equations to simple algebra system and solve the system.
"""
struct FODEMatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
@article{2009,
   title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
   DOI={10.1016/j.jcp.2009.01.014},
   author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}
"""
"""
Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional partial differential equations to simple algebra system and solve the system.
"""
struct FPDEMatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
@inproceedings{Podlubny1998FractionalDE,
  title={Fractional differential equations},
  author={Igor Podlubny},
  year={1998}
}
"""



"""
    solve(equation, right, h, T, MatrixDiscrete())

Using the **Matrix Discretization algorithm** proposed by [Prof Igor Podlubny](http://people.tuke.sk/igor.podlubny/index.html) to approximate the numerical solution.
"""
function solve(equation, right, h, T, ::FODEMatrixDiscrete)
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


"""

    solve(α, β, T, M, N, ::FPDEProblem)

When using the Martix 


"""
function solve(α, β, T, M, N, ::FPDEMatrixDiscrete)
    h=T/(M-1)
    τ=h^2/6
    TMatrix = kron(D(N-1, α, τ)', zeros(M, M)+I)
    SMatrix = kron(zeros(N-1, N-1)+I, RieszMatrix(β, M, h))

    system = TMatrix-SMatrix

    # Handling boundary conditions
    BMatrix = kron(zeros(N-1, N-1)+I, eliminator(M, [1 M]))
    system = system*BMatrix'

    left = BMatrix*system

    result = left\ones(size(left, 1), 1)
    return result

end
function RieszMatrix(α, N, h)
    caputo=B(N+1, α)
    caputo=caputo[2:(N+1), 1:N]
    result=1/2*(caputo+caputo')
    result=h^(-α)*result

    return result
end
function B(N, p)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i in range(1, N, step=1)
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end


"""
    diffusion(α, β)

**Diffusion equation** with time fractional derivative.
"""
function diffusion(α, β)
    solve(α, β, 1, 21, 148, FPDEMatrixDiscrete())
end