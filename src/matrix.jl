import FractionalDiffEq: solve, FractionalDiffEqAlgorithm

using LinearAlgebra, InvertedIndices

"""
Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional ordinary differential equations to simple algebra system and solve the system.

@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
"""
struct FODEMatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional partial differential equations to simple algebra system and solve the system.

@article{2009,
   title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
   DOI={10.1016/j.jcp.2009.01.014},
   author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}
"""
struct FPDEMatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
    solve(equation, right, h, T, MatrixDiscrete())

Using the **Matrix Discretization algorithm** proposed by [Prof Igor Podlubny](http://people.tuke.sk/igor.podlubny/index.html) to obtain the numerical solution.

## References

@inproceedings{Podlubny1998FractionalDE,
  title={Fractional differential equations},
  author={Igor Podlubny},
  year={1998}
}
"""
function solve(equation, right, highestorder, h, T, ::FODEMatrixDiscrete)
    
    
    N=Int64(T/h)
    rows = collect(1:highestorder)
    equation = eliminator(N, rows)*equation*eliminator(N, rows)'

    if typeof(right) <: Number
        rightside = eliminator(N, rows)*right*ones(N)
    else
        rightside = eliminator(N, rows)*right.(collect(h:h:T))
    end

    result = equation\rightside
    return vcat(zeros(highestorder), result)
end

"""
    eliminator(n, row)

Compute the eliminator matrix Sₖ by omiting n-th row
"""
function eliminator(n, row)
    temp = zeros(n, n) + I
    return temp[Not(row), :]
end

"""
Generating elements in Triangular Strip Matrix.
"""
function omega(n, α)
    omega = zeros(n+1)

    omega[1]=1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1] = (1 - (α+1)/i)*omega[i]
    end
    
    return omega

end

"""
    D(N, α, h)

Using D function to construct left hand equations.

### Example

Suppose we have a equation:

```math
y''(t)+D^{frac{3}{2}}_t y(t)=1
```

To represent the left side equation, use ``D(size, order, step)`` to construct the derivative part.

```julia-repl
equation=D(100, 2, 0.01)+D(100, 3/2, 0.01)
```

Then the right side can be constructed like:

```julia-repl
right = ones(N)
```

After we construct both sides, all we need to do is:

```julia-repl
solve(equation, right, h, T)
```

And then we can get the numerical approximation.

!!! info
    Here ```N``` is the size of discrete matrix.
"""
function D(N, α, h)
    result = zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i in range(1, N, step=1)
        result[i, 1:i] = reverse(temp[1:i])
    end

    return h^(-α)*result
end

function F(N, α, h)
    result=zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i ∈ 1:N
        result[i, 1:i]=reverse(temp[1:i])
    end

    result=reverse(reverse(result, dims=1), dims=2)

    return (-1)^(ceil(α))*h^(-α)*result
end

"""
    bagleytorvik(p1, p2, p3, T, h)

By specifying the parameters of Bagley Torvik Equation, we can use **bagleytorvik** to directly obtain the numerical approximation of a bagley torvik equation.

!!! info "p2 ≠ 0"
    Please note that the parameter of fractional derivative part must not be 0
"""
function bagleytorvik(p1, p2, p3, right, T, h)
    N=Int64(T/h+1)
    equation = p1*D(N, 2, h)+p2*D(N, 1.5, h)+p3*(zeros(N, N)+I)
    equation = eliminator(N, [1,2])*equation*eliminator(N, [1,2])'
    
    if typeof(right) <: Number
        rightside = eliminator(N, rows)*right*ones(N)
    else
        rightside = eliminator(N, rows)*right.(collect(h:h:T))
    end
    
    result = equation\rightside

    return vcat(zeros(2), result)
end


"""

    solve(α, β, T, M, N, FPDEMatrixDiscrete())

When using the Martix 
"""
function testsolve(α, β, T, M, N)
    h = T/(M-1)
    τ = h^2/6
    TMatrix = kron(D(N-1, α, τ)', zeros(M, M) + I)
    SMatrix = kron(zeros(N-1, N-1) + I, RieszMatrix(β, M, h))

    system = TMatrix-SMatrix

    # Handling boundary conditions
    BMatrix = kron(zeros(N-1, N-1) + I, eliminator(M, [1 M]))
    system = system*BMatrix'

    left = BMatrix*system

    result = left\(8*ones(size(left, 1), 1))

    return result
end

#=

tmp = testsolve(0.7, 0.8, 1, 21, 148)

YS = reshape(tmp, 19, 147)
Y = reverse(YS, dims=2)
U = 5 .*copy(Y)


rows, columns = size(U)

temp = [zeros(1, columns); U; zeros(1, columns)]
U=[zeros(1, 21)' temp]

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)
print(size(XX))

using Plots

plot3d(XX, YY, U)

=#

function meshgrid(xin,yin)
    nx=length(xin)
    ny=length(yin)
    xout=zeros(ny,nx)
    yout=zeros(ny,nx)
    for jx=1:nx
        for ix=1:ny
            xout[ix,jx]=xin[jx]
            yout[ix,jx]=yin[ix]
        end
    end
    return (x=xout, y=yout)
end


# Construct Riesz Symmetric Matrix
function RieszMatrix(α, N, h)
    caputo = B(N+1, α)
    caputo = caputo[2:(N+1), 1:N]
    result = 1/2*(caputo+caputo')
    result = h^(-α)*result

    return result
end
function B(N, p)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end


"""
    diffusion(α, β)

**Diffusion equation** with time fractional derivative.

α and β are the coefficients of the time derivative and spatial derivative.
"""
function diffusion(α, β)
    solve(α, β, 1, 21, 148, FPDEMatrixDiscrete())
end



## An fractional partial differential equation Example



