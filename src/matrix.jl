import FractionalDiffEq.solve

using LinearAlgebra, InvertedIndices, Plots


abstract type FractionalDiffEqAlgorithm end


struct MatrixDiscrete <: FractionalDiffEqAlgorithm end


"""
@inproceedings{Podlubny2000MATRIXAT,
  title={MATRIX APPROACH TO DISCRETE FRACTIONAL CALCULUS},
  author={Igor Podlubny},
  year={2000}
}
"""


"""
    eliminator(n, row)

Compute the eliminator matrix Sₖ by omiting n-th row
"""
function eliminator(n, row)
    temp = zeros(n, n)+I
    return temp[Not(row), :]
end



"""
@inproceedings{Podlubny1998FractionalDE,
  title={Fractional differential equations},
  author={Igor Podlubny},
  year={1998}
}
"""

"""
"""
function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    for i in range(1, n, step=1)
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega

end

function B(N, p, h)
    result=zeros(N, N)
    temp=omega(N, p)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    return h^(-p)*result
end

function F(N, p, h)
    result=zeros(N, N)
    temp = omega(N, p)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    result=reverse(reverse(result, dims=1), dims=2)

    return (-1)^(ceil(p))*h^(-p)*result
end


function solve(p1, α, p2, c, h, T, ::MatrixDiscrete)
    
    n=Int64(floor(α))
    rows=collect(1:n)
    N=Int64(T/h+1)
    equation = p1*B(N, α, h)+p2*(zeros(N, N)+I)
    equation = eliminator(N, rows)*equation*eliminator(N, rows)'
    righthand = c*eliminator(N, rows)*ones(N)
    result = equation\righthand
    return vcat(zeros(n), result)
end