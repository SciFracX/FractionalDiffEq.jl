function Rmat(n, l)
    if l == 0
        e = ones(n)
        R = myspdiagm(spdiagm(-1 => 0.5*e, 0 => e, 1 => 0.5*e), n)
        R[2, 1] = 1
    else
        nn = collect(0:n-1)
        e = ones(n)
        e1 = 0.5*(nn.+1)./(nn.+l)
        e2 = 0.5*(nn.+2*l.-1)./(nn.+l)
        R = myspdiagm(spdiagm(-1 => e1, 0 => e, 1 => e2[2:end]), n)
    end
    return R
end

myspdiagm(M, N) = size(M, 1) !== N ? (return M[1:N, 1:N]) : (return M)