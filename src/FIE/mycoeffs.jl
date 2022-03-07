function mycoeffs(f, n, lam)
    F = Fun(f, Ultraspherical(lam))
    m = coefficients(F)
    if length(m) < n
        return [m; zeros(n-length(m))]
    else
        return m[1:n]
    end
end