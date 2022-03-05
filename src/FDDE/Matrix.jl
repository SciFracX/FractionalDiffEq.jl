"""
# Usage

    solve(limit, alpha, A, B, f, t0, x0, T, tau, h, MatrixForm())

### Reference

https://github.com/mandresik/system-of-linear-fractional-differential-delayed-equations
"""
struct MatrixForm <: FractionalDiffEqAlgorithm end

function solve(limit, alpha, A, B, f, t0, x0, T, tau, h, ::MatrixForm)
    var_num = length(A[:, 1])
    m = ceil(Int, alpha)

    t_tau = t0-tau
    temp = collect(t0-h:-h:t_tau)
    t = [temp[end:-1:1]; collect(t0:h:T)]

    index = length(temp)
    N = length(t)-index-1

    x = zeros(var_num, length(t))
    b = zeros(1, N)
    F = zeros(var_num, N)
    x0TP = x0(t0)

    @fastmath @inbounds @simd for i=1:index+1
        x[:, i] = x0(t[i])*[1; zeros(m-1, 1)]
    end



    function function_values(M, limit)
        cols::Int64 = size(M, 2)
        risk_index = zeros(var_num)
        col = 0
        M0 = copy(M)

        @fastmath @inbounds @simd for r=1:var_num
            @fastmath @inbounds @simd for s=1:cols
                if M0[r, s]==Inf || M0[r, s]==-Inf
                    risk_index = hcat(risk_index[:, 1:col], [r; s])
                    col += 1
                end
            end
        end
        Mt = zeros(var_num, cols*length(t))

        for j=1:length(t)
            Mt[:, (cols*j-cols+1):(cols*j)] = M
        end

        if risk_index != zeros(var_num, 1)
            risk_start_index = judgeeqsum(t, -1)
            inv_num = judgesum(t, 1)-judgeeqsum(t, -1)
            risks_num = size(risk_index, 2)
            @fastmath @inbounds @simd for j=1:inv_num
                @fastmath @inbounds @simd for k=1:risks_num
                    id = risk_start_index*cols + risk_index[2, k] + cols*(j-1)
                    if Mt[risk_index[1, k], id] > limit
                        Mt[risk_index[1, k], id] = limit
                    elseif Mt[risk_index[1, k], id] < -limit
                        M[risk_index[1, k], id] = -limit
                    end
                end
            end
        end
        return Mt
    end

    At = function_values(A, limit)
    Bt = function_values(B, limit)
    ft = function_values(f, limit)

    @fastmath @inbounds @simd for n=1:N
        col_num = size(A, 2)
        ind = col_num*(n+index)
        b[n] = (n+1)^alpha-n^alpha
        F[:, n] = At[:, ind-col_num+1:ind]*x[:, n+index] + Bt[:, ind-col_num+1:ind]*x[:, n]+ft[:, n+index]
        TP = zeros(var_num, 1)
        @fastmath @inbounds @simd for k=0:m-1
            TP += ((t[n+index+1]-t0)^k).*x0TP[:, k+1]./factorial(k)
        end
        sum = zeros(var_num, 1)
        @fastmath @inbounds @simd for j=1:n
            sum += b[n-j+1].*F[:, j]
        end
        x[:, n+index+1] = TP + h^alpha.*sum
    end
    return x        
end

function judgesum(t, thres)
    result = 0
    for i in t
        i < thres ? result+=i : continue
    end
    return result
end

function judgeeqsum(t, thres)
    result = 0
    for i in t
        i <= thres ? result+=i : continue
    end
    return result
end