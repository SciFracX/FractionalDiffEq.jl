using SparseArrays

function Qmat(N, m, lam)
    # Integration matrix
    Q = 1
    if round(m) !== m
        Q = Qmat05(N, lam)     
    end
    Q1 = Qmat1(N, lam)
    Q = Q*Q1
    return Q
end
    
function Qmat05(N, lam)
    # Half integral
    if lam == 0.5
        v = (2/sqrt(pi)) ./ (2*collect(0:N-1).+1)
        Q = myspdiagm(spdiagm(0=>v, 1=>-v[2:end]), N)
        #Q = spdiags([v -v], [0 1], N, N)
    elseif lam == 1
        v = (sqrt(pi)/2) .* ones(N)
        Q = myspdiagm(spdiagm(-1=>v[2:end], 0=>v), N)
        #Q = spdiags([v, v], [-1, 0], N, N)
    end
    return Q
end
    
function Qmat1(N, lam)
    # Full integral
    if lam == 0.5
        v = 1 ./ (2*collect(0:N-1).+1)
        Q = myspdiagm(spdiagm(-1 => v, 1 => -v[2:end]), N)
        #Q = spdiags([v,-v], [-1, 1], N, N)
        Q[1, 1] = 1
    elseif lam == 1
        nn = collect(0:N)
        v = 1 ./ (2*nn.+1)
        v2 = 2 ./ (4*nn.^2 .-1)
        Q = myspdiagm(spdiagm(-1=>v[2:end], 0=>v2[2:end], 1=>-v[2:end-1]), N)
        #Q = spdiags([v[2:end], v2[2:end], -v[1:end-1]], [-1:1], N, N);
    end
    return Q
end

myspdiagm(M, N) = size(M, 1) !== N ? (return M[1:N, 1:N]) : (return M)