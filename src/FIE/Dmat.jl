function Dmat(N, m, lam)

    # Initialise D with identity or half-derivative accordingly:
    if (round(m) == m) && (lam == 1/2)
        e = (2^m*gamma(m+1/2)/sqrt(pi)) .* ones(N)
        D = myspdiagm(spdiagm(m => e), N)
    elseif (round(m) !== m) && (lam == 1)
        e = 2^(m-1/2)*gamma(m+1) * ones(N,1)
        D = myspdiagm(spdiagm(Int(m-1/2) => e, Int(m+1/2) =>e), N)
    elseif round(m) == m
        D = 1
        for k = 0:floor(m)-1
            D = Dmat1(N, -k+0.5, lam+k) * D
        end
    else
        D = Dmat05(N, lam)
        for k = 0:floor(m)-1
            D = Dmat1(N, -k-0.5, lam+k+0.5) * D
        end
    end
    return D
end
    
function Dmat05(N, lam)
    # Half derivative: d^{1/2}/dx^{1/2} (1+x)^(lam-1/2)*C^{(lam)}_n(x)
    e = (gamma(lam+0.5)/gamma(lam)) .* ones(N)
    D = myspdiagm(spdiagm(0 => e, 1 => e), N)
    return D
end

function Dmat1(N, mu, lam)
    # Full derivative: d/dx (1+x)^mu*C^{(l)}_n(x)
    e = ones(N)
    if mu == 0
        D = 2*lam*myspdiagm(spdiagm(1 => e), N)
    else
        v = (mu-lam)./(collect(0:N-1).+lam)
        D = myspdiagm(spdiagm(0 => lam.*(1 .+v), 1 => lam.*2 .*e, 2 => (lam.*(1 .-v))[3:end]), N)
    end
end

myspdiagm(M, N) = size(M, 1) !== N ? (return M[1:N, 1:N]) : (return M)