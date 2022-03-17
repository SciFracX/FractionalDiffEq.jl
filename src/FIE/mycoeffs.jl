# Generate right hand coefficients
function mycoeffs(f, n, lam)
    F = Fun(f, Ultraspherical(lam))
    m = coefficients(F)
    length(m) < n ? (return [m; zeros(n-length(m))]) : (return m[1:n])
end