function mycoeffs(f, n, lam)
    F = Fun(f, Ultraspherical(lam))
    return coefficients(F)[1:n]
end