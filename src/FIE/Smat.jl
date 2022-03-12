function Smat(n, lam)
    if lam == 0
        e = ones(n)
        S = myspdiagm(spdiagm(0 => 0.5*e, 2 => -0.5*e), n)
        S[1] = 1
    else
        e = lam./(collect(0:n-1).+lam)
        S = myspdiagm(spdiagm(0 => e, 2 => -e), n)
    end
    return S
end

myspdiagm(M, N) = size(M, 1) !== N ? (return M[1:N, 1:N]) : (return M)