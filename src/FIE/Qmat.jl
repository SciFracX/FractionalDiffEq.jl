function Qmat(N, m, lam)
    # Integration matrix
    Q = 1
    round(m) !== m ? (return Q = Qmat05(N, lam)) : nothing
    Q = Q*Qmat1(N, lam)
    return Q
end

function Qmat05(N, lam)
    # Half integral
    if lam == 0.5
        v = (2/sqrt(pi)) ./ (2*collect(0:N-1).+1)
        Q = myspdiagm(spdiagm(0=>v, 1=>-v[2:end]), N)
    elseif lam == 1
        v = (sqrt(pi)/2) .* ones(N)
        Q = myspdiagm(spdiagm(-1=>v[2:end], 0=>v), N)
    end
    return Q
end
    
function Qmat1(N, lam)
    # Full integral
    if lam == 0.5
        v = 1 ./ (2*collect(0:N-1).+1)
        Q = myspdiagm(spdiagm(-1 => v, 1 => -v[2:end]), N)
        Q[1, 1] = 1
    elseif lam == 1
        nn = collect(0:N)
        v = 1 ./ (2*nn.+1)
        v2 = 2 ./ (4*nn.^2 .-1)
        Q = myspdiagm(spdiagm(-1=>v[2:end], 0=>v2[2:end], 1=>-v[2:end-1]), N)
    end
    return Q
end

myspdiagm(M, N) = size(M, 1) !== N ? (return M[1:N, 1:N]) : (return M)