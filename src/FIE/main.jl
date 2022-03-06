function solve(f, e, n)
    # Initialize operators
    Z = spzeros(n, n)
    I1 = sparse(Matrix(I, n, n))
    II = [I1 Z; Z I1]
    Q05 = lam -> Qmat(n, 0.5, lam)
    QQ05 = [Z Q05(1); Q05(0.5) Z]
    QQ = m->QQ05^(2*m)

    rhs = [mycoeffs(e, n, 0.5); mycoeffs(f, n, 1)]
    idx = [collect(1:n); collect(n+1:2*n)]
    A = A[idx, idx]

    # Initialize solution
    u=zeros(size(A, 1))
    u[idx] = A\rhs
end