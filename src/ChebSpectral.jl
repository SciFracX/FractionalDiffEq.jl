"""
    solve(prob::SingleTermFODEProblem, n, ChebSpectral())

Using Chebyshev spectral method to solve single term fractional ordinary differential equations.


!!! warning
    Chenyshev spectral method only support for linear fractional differential equations, and the time span of the problem should be ``[-1, 1]``.

### References

```tex
@thesis{
    title = {Spectral Methods for Fractional Differential Equations}
    author = {Chinenye Assumpta Nnakenyi, J. A. C. Weideman, Nicholas Hale}
    year = {2015}
}
```
"""
struct ChebSpectral <: FractionalDiffEqAlgorithm end
#FIXME: Fix the initial conditions

function solve(FODE::SingleTermFODEProblem, n, ::ChebSpectral)
    @unpack f, α, u0, T = FODE
    (D, x) = Cheb(n, α)
    result = D\f.(x)
    return result
end

function VandermondeLegendre(n)
    x = [-cos(π*j/n) for j in range(0, n, step=1)]
    A = x*x'
    A[1, :] = ones(n+1)
    A[2, :] = x
    for k = 1:n-1
        for i = 0:n
            A[k+2, i+1] = ((2*k+1)*A[2, i+1]*A[k+1, i+1] - k*A[k, i+1])/(k+1)
        end
    end
    return A'
end

function VandermondeJacobi(n, mu)
    x = [-cos(π*j/n) for j in range(0, n, step=1)]
    A = x*x'
    a = -mu
    b = mu
    A0 = 0.5*(a+b)+1
    B0 = 0.5*(a-b)
    A[1, :] = ones(n+1)
    A[2, :] = A0.*x .+ B0
    for k = 1:n-1
        for i = 0:n
            An =((2*k+a+b+1)*(2*k+a+b+2))/(2*(k+1)*(k+a+b+1))
            Bn = ((a^2-b^2)*(2*k+a+b+1))/(2*(k+1)*(k+a+b+1)*(2*k+a+b))
            Cn = ((k+a)*(k+b)*(2*k+a+b+2))/((k+1)*(k+a+b+1)*(2*k+a+b))
            A[k+2, i+1] = (An*x[i+1]+Bn)*A[k+1, i+1] - Cn*A[k, i+1]
        end
    end
    return A'
end

"""
Chebyshev fractional spectral differentiation matrix
"""
function Cheb(n, mu)
    if mu == 1
        if n == 0
            @assert("Wrong")
        else
            k = LinRange(0, 1, n+1)
            D = k*k'
            D[1, 1]= -(2*n^2 +1)/6
            D[n, n]=(2*n^2 +1)/6
            D[n, 1] = (-1)^n/2
            D[1, n] = -(-1)^n/2
            for j = 0:n
                for l = 0:n
                    if j !== l
                        if k[j+1] == 0 || k[j+1]==n
                            c=2
                        else
                            c = 1
                            D[j+1, l+1] = -c*(-1)^(j+l)/(c*(cos(k[j+1]*π) - cos(k[l+1]*π)))
                        end
                    end
                end
            end

            for i = 1:n-1
                D[i+1, i+1] = -(-1)*np.cos(k[i]*pi)/(2*(1-(cos(k[i+1]*pi))^2))
                D[1, i+1] = -2*(-1)^(i)/(1 - (cos(k[i+1]*pi)))
                D[n, i+1] = 2*(-1)^(i+n)/(1 + (cos(k[i+1]*pi)))
                D[i+1, 1] = 1*(-1)^(i)/(2*(1 - (cos(k[i+1]*pi))))
                D[i+1, n] = -(-1)^(i+n)/(2*(1 + (cos(k[i+1]*pi))))
            end
            x = [-cos(π*j/n) for j in range(0, n, step=1)]
            return D, x
        end
        else
            k = zeros(n)
            A = k*k'
            A1 = k*k'
            B = []
            x = [-cos(π*j/(n-1)) for j in range(0, n-1, step=1)]
            B1 = @. 1/((x+1)^mu)
            for k=1:n
                push!(B, gamma(k+mu)/gamma(k))
            end
            for i=0:n-1
                for j=0:n-1
                    if i == j
                        A[i+1, j+1] = B[i+1]
                        A1[i+1, j+1] = B1[i+1]
                    end
                end
            end
            P = VandermondeLegendre(n-1)
            j1 = inv(VandermondeJacobi(n-1, mu))
            D = P*A*j1*A1
            D = D[2:n, 2:n]
            x = x[2:end]
            return D, x
        end
end
#=
n=30
#tspan = collect(-1:0.01:1)

right(t)=t*exp(-t)
#RHS = rightfun.(tspan)
(D, x) = Cheb(n, 0.5)
sol = D\rightfun.(x)

using Plots
plot(x, sol)
=#
#=
U = @. (1+x)^(6+9/17)
plot!(x, U, ls=:dash)
=#