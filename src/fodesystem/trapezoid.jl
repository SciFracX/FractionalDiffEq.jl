"""
    solve(prob::FODEProblem, FLMMTrap())

Use [Trapezoidal](https://en.wikipedia.org/wiki/Trapezoidal_rule_(differential_equations)) with generating function ``f(x)=\\frac{1+x}{2(1-x)^\\alpha}`` generated weights fractional linear multiple steps method to solve system of FODE.

### References

```tex
@article{Garrappa2015TrapezoidalMF,
  title={Trapezoidal methods for fractional differential equations: Theoretical and computational aspects},
  author={Roberto Garrappa},
  journal={ArXiv},
  year={2015},
  volume={abs/1912.09878}
}
```
"""
struct FLMMTrap <: FODESystemAlgorithm end

function solve(prob::FODEProblem, h, ::FLMMTrap; reltol=1e-6, abstol=1e-6)
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    u0 = u0
    alpha = order[1]
    maxiters = 100

    # issue [#64](https://github.com/SciFracX/FractionalDiffEq.jl/issues/64)
    max_order = findmax(order)[1]
    if max_order > 1
        @error "This method doesn't support high order FDEs"
    end

    Jfdefun(t, u) = jacobian_of_fdefun(f, t, u, p)

    m_alpha::Int = ceil.(Int, alpha)
    m_alpha_factorial = factorial.(collect(0:m_alpha-1))
    # Structure for storing information on the problem
    
    problem_size = size(u0, 1)
    
    
    # Check compatibility size of the problem with size of the vector field
    
    # Number of points in which to evaluate the solution or the weights
    r::Int = 16
    N::Int = ceil(Int, (tfinal-t0)/h)
    Nr::Int = ceil(Int, (N+1)/r)*r
    Q::Int = ceil(Int, log2((Nr)/r))-1
    global NNr = 2^(Q+1)*r

    # Preallocation of some variables
    y = zeros(problem_size, N+1)
    fy = zeros(problem_size, N+1)
    zn = zeros(problem_size, NNr+1)

    # Evaluation of convolution and starting weights of the FLMM
    (omega, w, s) = TrapWeights(alpha, NNr+1)
    halpha = h^alpha
    
    # Initializing solution and proces of computation
    t = collect(0:N)*h
    y[:, 1] = u0[:, 1]
    temp = zeros(problem_size)
    f(temp, u0[:, 1], p, t0)
    fy[:, 1] = temp
    (y, fy) = TrapFirstApproximations(t, y, fy, abstol, maxiters, s, halpha, omega, w, problem_size, f, Jfdefun, u0, m_alpha, t0, m_alpha_factorial, p)
    (y, fy) = TrapTriangolo(s+1, r-1, 0, t, y, fy, zn, N, abstol, maxiters, s, w, omega, halpha, problem_size, f, Jfdefun, u0, m_alpha, t0, m_alpha_factorial, p)
    
    # Main process of computation by means of the FFT algorithm
    nx0 = 0; ny0 = 0
    ff = zeros(1, 2^(Q+2), 1)
    ff[1:2] = [0 2]
    for q = 0:Q
        L::Int = 2^q
        (y, fy) = TrapDisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, zn, N, abstol, maxiters, s, w, omega, halpha, problem_size, f, Jfdefun, u0, m_alpha, t0, m_alpha_factorial, p)
        ff[1:4*L] = [ff[1:2*L]; ff[1:2*L-1]; 4*L]
    end
    # Evaluation solution in TFINAL when TFINAL is not in the mesh
    if tfinal < t[N+1]
        c = (tfinal - t[N])/h
        t[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end
    t = t[1:N+1]; y = y[:, 1:N+1]
    return FODESystemSolution(t, y)
end


function TrapDisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, zn, N::Int, abstol, maxiters, s, w, omega, halpha, problem_size, fdefun, Jfdefun, y0, m_alpha, t0, m_alpha_factorial, p)
    nxi::Int = copy(nx0); nxf::Int = copy(nx0 + L*r - 1)
    nyi::Int = copy(ny0); nyf::Int = copy(ny0 + L*r - 1)
    is::Int = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[is] = nxi; s_nxf[is] = nxf ; s_nyi[is] = nyi; s_nyf[is] = nyf
    i_triangolo::Int = 0;  stop = false
    while ~stop
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1)
        zn = TrapQuadrato(nxi, nxf, nyi, nyf, fy, zn, omega, problem_size)
        (y, fy) = TrapTriangolo(nxi, nxi+r-1, nxi, t, y, fy, zn, N, abstol, maxiters, s, w, omega, halpha, problem_size, fdefun, Jfdefun, y0, m_alpha, t0, m_alpha_factorial, p)
        i_triangolo = i_triangolo + 1
        
        if ~stop
            if nxi+r-1 == nxf   # Il triangolo finisce dove finisce il quadrato -> si scende di livello
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is] + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            else # Il triangolo finisce prima del quadrato -> si fa un quadrato accanto
                nxi = nxi + r; nxf = nxi+r-1; nyi = nyf+1; nyf = nyf+r
                is = is + 1
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            end
        end
        
    end
    return y, fy
end

function TrapQuadrato(nxi, nxf, nyi, nyf, fy, zn, omega, problem_size)
    coef_beg::Int = nxi-nyf; coef_end::Int = nxf-nyi+1
    funz_beg::Int = nyi+1; funz_end::Int = nyf+1
    vett_coef = omega[coef_beg+1:coef_end+1]
    vett_funz = [fy[:, funz_beg:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
    zzn = real(fast_conv(vett_coef, vett_funz))
    zn[:, nxi+1:nxf+1] = zn[:, nxi+1:nxf+1] + zzn[:, nxf-nyf:end-1]
    return zn
end

function TrapTriangolo(nxi, nxf, j0, t, y, fy, zn, N, abstol, maxiters, s, w, omega, halpha, problem_size, fdefun, Jfdefun, y0, m_alpha, t0, m_alpha_factorial, p)
    for n = nxi:min(N, nxf)
        n1::Int = n+1
        St = TrapStartingTerm(t[n1], y0, m_alpha, t0, m_alpha_factorial)
        
        Phi = zeros(problem_size, 1)
        for j = 0:s
            Phi = Phi + w[j+1, n1]*fy[:, j+1]
        end
        for j = j0:n-1
            Phi = Phi + omega[n-j+1]*fy[:, j+1]
        end
        Phi_n = St + halpha*(zn[:, n1] + Phi)
        
        yn0 = y[:, n]
        temp = zeros(length(yn0))
        fdefun(temp, yn0, p, t[n1])
        fn0 = temp#f_vectorfield(t[n1], yn0, fdefun)
        Jfn0 = Jf_vectorfield(t[n1], yn0, Jfdefun)
        Gn0 = yn0 - halpha*omega[1]*fn0 - Phi_n
        stop = false; it::Int = 0
        while ~stop            
            JGn0 = zeros(problem_size, problem_size)+I - halpha*omega[1]*Jfn0
            global yn1 = yn0 - JGn0\Gn0
            global fn1 = zeros(length(yn1))#f_vectorfield(t[n1], yn1, fdefun)
            fdefun(fn1, yn1, p, t[n1])
            Gn1 = yn1 - halpha*omega[1]*fn1 - Phi_n
            it = it + 1
            
            stop = (norm(yn1-yn0, Inf) < abstol) || (norm(Gn1, Inf)<abstol)
            if it > maxiters && ~stop
                @warn "Non Convergence"
                stop = true
            end
            
            yn0 = yn1; Gn0 = Gn1
            if ~stop
                Jfn0 = Jf_vectorfield(t[n1], yn0, Jfdefun)
            end
            
        end
        y[:, n1] = yn1
        fy[:, n1] = fn1
    end
    return y, fy
end

function TrapFirstApproximations(t, y, fy, abstol, maxiters, s, halpha, omega, w, problem_size, fdefun, Jfdefun, y0, m_alpha, t0, m_alpha_factorial, p)
    m = problem_size
    Im = zeros(m, m)+I ; Ims = zeros(m*s, m*s)+I
    Y0 = zeros(s*m, 1); F0 = copy(Y0); B0 = copy(Y0)
    for j = 1 : s
        Y0[(j-1)*m+1:j*m, 1] = y[:, 1]
        temp = zeros(length(y[:, 1]))
        fdefun(temp, y[:, 1], p, t[j+1])
        F0[(j-1)*m+1:j*m, 1] = temp#f_vectorfield(t[j+1], y[:, 1], fdefun)
        St = TrapStartingTerm(t[j+1], y0, m_alpha, t0, m_alpha_factorial)
        B0[(j-1)*m+1:j*m, 1] = St + halpha*(omega[j+1]+w[1, j+1])*fy[:, 1]
    end
    W = zeros(s, s)
    for i = 1:s
        for j = 1:s
            if i >= j
                W[i, j] = omega[i-j+1] + w[j+1, i+1]
            else
                W[i, j] = w[j+1, i+1]
            end
        end
    end
    W = halpha*kron(W, Im)
    G0 = Y0 - B0 - W*F0
    JF = zeros(s*m, s*m)
    for j = 1:s
        JF[(j-1)*m+1:j*m, (j-1)*m+1:j*m] = Jf_vectorfield(t[j+1], y[:, 1], Jfdefun)
    end
    stop = false; it::Int = 0
    F1 = zeros(s*m, 1)
    while ~stop
        JG = Ims - W*JF
        global Y1 = Y0 - JG\G0
        
        for j in 1:s
            temp = zeros(length(Y1[(j-1)*m+1:j*m, 1]))
            fdefun(temp, Y1[(j-1)*m+1:j*m, 1], p, t[j+1])
            F1[(j-1)*m+1:j*m, 1] = temp#f_vectorfield(t[j+1], Y1[(j-1)*m+1:j*m, 1], fdefun)
        end
        G1 = Y1 - B0 - W*F1
        
        it = it + 1
        
        stop = (norm(Y1-Y0, Inf) < abstol) || (norm(G1, Inf) <  abstol)
        if it > maxiters && ~stop
            @warn "Non Convergence"
            stop = 1
        end
        
        Y0 = Y1 ; G0 = G1
        if ~stop
            for j = 1 : s
                JF[(j-1)*m+1:j*m, (j-1)*m+1:j*m] = Jf_vectorfield(t[j+1], Y1[(j-1)*m+1:j*m, 1], Jfdefun)
            end
        end
        
    end
    for j = 1 : s
        y[:, j+1] = Y1[(j-1)*m+1:j*m, 1]
        fy[:, j+1] = F1[(j-1)*m+1:j*m, 1]
    end
    return y, fy
end

function TrapWeights(alpha, N)
    # Trapezoidal method with generating function ((1+x)/2/(1-x))^alpha
    omega1 = zeros(1, N+1); omega2 = copy(omega1)
    omega1[1] = 1; omega2[1] = 1
    alpha_minus_1 = alpha - 1 ; alpha_plus_1 = alpha + 1
    for n = 1 : N
        omega1[n+1] = (alpha_plus_1/n - 1)*omega1[n]
        omega2[n+1] = (1 + alpha_minus_1/n)*omega2[n]
    end
    x = fft([omega1 zeros(size(omega1))])
    y = fft([omega2 zeros(size(omega2))])
    omega = ifft(x.*y) 
    omega = omega[1:N+1]/2^alpha
    omega = real.(omega)

    k = floor(1/abs(alpha))
    if abs(k - 1/alpha) < 1.0e-12
        A = collect(0:k)*abs(alpha)
    else
        A = [collect(0:k)*abs(alpha); 1]
    end
    s = length(A) - 1
    # Generation of the matrix and the right hand--side vectors of the system
    nn = collect(0:N)
    V = zeros(s+1, s+1); jj_nu = zeros(s+1, N+1); nn_nu_alpha = copy(jj_nu)
    for i = 0:s
        nu = A[i+1]
        V[i+1, :] = collect(0:s).^nu
        jj_nu[i+1, :] = nn.^nu
        if alpha > 0
            nn_nu_alpha[i+1, :] = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha)
        else
            if i == 0
                nn_nu_alpha[i+1, :] = zeros(1, N+1)
            else
                nn_nu_alpha[i+1, :] = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha)
            end
        end
    end

    temp = fast_conv([omega  zeros(size(omega))], [jj_nu zeros(size(jj_nu))])
    temp = real.(temp)
    b = nn_nu_alpha - temp[:, 1:N+1]
    # Solution of the linear system with multiple right-hand side
    w = real.(V\b)

    return omega, w, s
end

function TrapStartingTerm(t,y0, m_alpha, t0, m_alpha_factorial)
    ys = zeros(size(y0, 1), 1)
    for k = 1:m_alpha
        ys = ys + (t-t0)^(k-1)/m_alpha_factorial[k]*y0[:, k]
    end
    return ys
end