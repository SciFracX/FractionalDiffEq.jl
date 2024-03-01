@concrete struct MatrixDiscreteCache{T}
    prob
    alg
    mesh

    highest_order

    equation
    rightside

    kwargs
end

Base.eltype(::MatrixDiscreteCache{T}) where {T} = T

function SciMLBase.__init(prob::MultiTermsFODEProblem, alg::MatrixDiscrete; dt = 0.0, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack parameters, orders, f, u0, tspan, p = prob
    T = eltype(u0)
    t0 = tspan[1]; tfinal = tspan[2]
    mesh = collect(T, t0+dt:dt:tfinal)
    N::Int64 = ceil(Int, (tfinal-t0)/dt)
    highest_order = findmax(ceil.(Int, orders))[1]
    rows = collect(Int64, 1:highest_order)

    equation = zeros(T, N, N)

    for (i, j) in zip(parameters, orders)
        equation += i*D(N, j, dt)
    end


    equation = eliminator(N, rows)*equation*eliminator(N, rows)'

    if typeof(f) <: Number # Handling right hand side
        rightside = eliminator(N, rows)*f*ones(T, N)
    else
        tmp = map(x -> f(nothing, p, x), mesh)
        rightside = eliminator(N, rows)*(tmp .+ ic_handling(orders, parameters, u0))
    end

    return MatrixDiscreteCache{T}(prob, alg, mesh, highest_order, equation, rightside, kwargs)
end

function SciMLBase.solve!(cache::MatrixDiscreteCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, highest_order, equation, rightside, kwargs = cache
    result = equation\rightside
    result = vcat(zeros(highest_order), result)
    y = collect(Vector{eltype(result)}, eachrow(result))

    return DiffEqBase.build_solution(prob, alg, mesh, y)
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

    omega[1] = 1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1] = (1 - (α+1)/i)*omega[i]
    end
    
    return omega
end

"""
    D(N, α, dt)

Using D function to construct left hand side equations.

!!! info
    Here ```N``` is the size of discrete matrix.
"""
function D(N, α, dt)
    α == 0 ? (return zeros(N, N) + I) : nothing # When α=0, D is an identity matrix.
    result = zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i in range(1, N, step=1)
        result[i, 1:i] = reverse(temp[1:i])
    end

    return dt^(-α)*result
end

function F(N, α, dt)
    result = zeros(N, N)
    temp = omega(N, α)

    @fastmath @inbounds @simd for i ∈ 1:N
        result[i, 1:i] = reverse(temp[1:i])
    end

    result = reverse(reverse(result, dims=1), dims=2)

    return (-1)^(ceil(α))*dt^(-α)*result
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
function RieszMatrix(α, N, dt)
    caputo = B(N+1, α)
    caputo = caputo[2:(N+1), 1:N]
    result = 1/2*(caputo+caputo')
    result = dt^(-α)*result

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
