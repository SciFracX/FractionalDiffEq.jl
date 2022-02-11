"""
# Usage
    solve(prob::MultiTermsFODEProblem, h, T, FODEMatrixDiscrete())

Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional ordinary differential equations to simple algebra system and solve the system.

### References

```tex
@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
```
"""
struct FODEMatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
    solve(α, β, T, M, N, FPDEMatrixDiscrete())

Using [triangular strip matrices](https://en.wikipedia.org/wiki/Triangular_matrix) to discrete fractional partial differential equations to simple algebra system and solve the system.

```tex
@article{2009,
   title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
   DOI={10.1016/j.jcp.2009.01.014},
   author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}
```
"""
struct FPDEMatrixDiscrete <: FractionalDiffEqAlgorithm end



function solve(prob::MultiTermsFODEProblem, h, T, ::FODEMatrixDiscrete)
    leftparameters, leftorders, right = prob.parameters, prob.orders, prob.rightfun
    N::Int64 = floor(Int, T/h)
    highestorder = Int64(findmax(ceil.(leftorders))[1])
    rows = collect(1:highestorder)

    equation = zeros(N, N)

    for (i, j) in zip(leftparameters, leftorders)
        equation += i*D(N, j, h)
    end


    equation = eliminator(N, rows)*equation*eliminator(N, rows)'

    if typeof(right) <: Number
        rightside = eliminator(N, rows)*right*ones(N)
    else
        rightside = eliminator(N, rows)*right.(collect(h:h:T))
    end

    result = equation\rightside
    result = vcat(zeros(highestorder), result)

    return result
end


"""
    eliminator(n, row)

Compute the eliminator matrix Sₖ by omiting n-th row
"""
function eliminator(n, row)
    temp = I(n)
    return temp[Not(row), :]
end

"""
Generating elements in Triangular Strip Matrix.
"""
function omega(n, α)
    omega = zeros(n+1)

    omega[1] = 1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1] = (1 - (α+1)/i)*omega[i]
    end
    
    return omega
end

"""
    D(N, α, h)

Using D function to construct left hand side equations.

!!! info
    Here ```N``` is the size of discrete matrix.
"""
function D(N, α, h)
    α == 0 ? (return zeros(N, N) + I) : nothing # When α=0, D is an identity matrix.
    result = zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i in range(1, N, step=1)
        result[i, 1:i] = reverse(temp[1:i])
    end

    return h^(-α)*result
end

function F(N, α, h)
    result = zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i ∈ 1:N
        result[i, 1:i] = reverse(temp[1:i])
    end

    result = reverse(reverse(result, dims=1), dims=2)

    return (-1)^(ceil(α))*h^(-α)*result
end

function meshgrid(xin,yin)
    nx = length(xin)
    ny = length(yin)
    xout = zeros(ny,nx)
    yout = zeros(ny,nx)
    for jx = 1:nx
        for ix = 1:ny
            xout[ix,jx] = xin[jx]
            yout[ix,jx] = yin[ix]
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
    result = zeros(N, N)
    temp = omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end

"""
    bagleytorvik(p1, p2, p3, T, h)

By specifying the parameters of Bagley Torvik Equation, we can use **bagleytorvik** to directly obtain the numerical approximation of a bagley torvik equation.

!!! info "p2 ≠ 0"
    Please note that the parameter of fractional derivative part must not be 0
"""
function bagleytorvik(p1, p2, p3, right, T, h)
    N = round(Int, T/h)+1
    equation = p1*D(N, 2, h) + p2*D(N, 1.5, h) + p3*(zeros(N, N)+I)
    equation = eliminator(N, [1,2])*equation*eliminator(N, [1,2])'

    rows = collect(1:2)
    
    if typeof(right) <: Number
        rightside = eliminator(N, rows)*right*ones(N)
    else
        rightside = eliminator(N, rows)*right.(collect(h:h:T))
    end
    
    result = equation\rightside

    return vcat(zeros(2), result)
end


function solve(α, β, T, M, N, ::FPDEMatrixDiscrete)
    h = T/(M-1)
    τ = h^2/6
    
    # Construct the time partial derivative matrix
    TMatrix = kron(D(N-1, α, τ)', zeros(M, M) + I)

    # Construct the spatial partial derivative matrix
    SMatrix = kron(zeros(N-1, N-1) + I, RieszMatrix(β, M, h))

    system = TMatrix-SMatrix

    # Handling boundary conditions
    BMatrix = kron(zeros(N-1, N-1) + I, eliminator(M, [1, M]))
    system = system*BMatrix'

    left = BMatrix*system

    tmp = left\(ones(size(left, 1), 1))

    YS = reshape(tmp, M-2, N-1)
    YS = reverse(YS, dims=2)
    U = YS
    rows, columns = size(U)

    # Boundry conditions
    U = [zeros(1, columns); U; zeros(1, columns)]
    result = [zeros(M) U]

    return result
end



#=
Test code for FPDEMatrixDiscrete, to explore the wide usage of FPDEMatrixDiscrete
=#
#=
function testsolve(α, β, T, M, N)
    h = T/(M-1)
    τ = h^2/6
    
    # Construct the time partial derivative matrix
    TMatrix = kron(D(N-1, α, τ)', zeros(M, M) + I)

    # Construct the spatial partial derivative matrix
    SMatrix = kron(zeros(N-1, N-1) + I, RieszMatrix(β, M, h))

    system = TMatrix+SMatrix

    #FIXME: Handling boundary conditions
    BMatrix = kron(zeros(N-1, N-1) + I, eliminator(M, [1, M]))
    system = system*BMatrix'

    leftside = BMatrix*system

    XX = range(1/2793, 1, step=1/2793)
    YY = range(1/2793, 1, step=1/2793)

    result = leftside\(-sin.(π .*XX).*sin.(π .*YY))


    return result
end

using Plots

using LinearAlgebra, InvertedIndices

tmp = testsolve(2, 2, 1, 21, 148)



YS = reshape(tmp, 19, 147)
YS = reverse(YS, dims=2)
U = YS


rows, columns = size(U)

U = [zeros(1, columns); U; zeros(1, columns)]
U=[zeros(1, 21)' U]

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)

plotlyjs()

plot(XX, YY, U, st=:surface)


=#


"""
    diffusion(α, β)

**Diffusion equation** with time fractional derivative.

α and β are the coefficients of the time derivative and spatial derivative.
"""
function diffusion(α, β)
    solve(α, β, 1, 21, 148, FPDEMatrixDiscrete())
end

## An fractional partial differential equation Example



