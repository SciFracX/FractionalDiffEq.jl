@concrete mutable struct BDFCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    halpha
    y
    fy
    zn
    p

    # FLMM coefficients
    m_alpha
    m_alpha_factorial
    r
    N
    Nr
    Q
    NNr

    # LMM weights
    omega
    w
    s

    dt
    reltol
    abstol
    maxiters
    
    kwargs
end

Base.eltype(::BDFCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::BDF;
                          dt=0.0, reltol=1e-6, abstol=1e-6, maxiters = 1000, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    all(x->x==order[1], order) ? nothing : throw(ArgumentError("BDF method is only for commensurate order FODE"))
    alpha = order[1] # commensurate ordre FODE
    (alpha > 1.0) && throw(ArgumentError("BDF method is only for order <= 1.0"))
    

    m_alpha = ceil.(Int, alpha)
    m_alpha_factorial = factorial.(collect(0:m_alpha-1))
    problem_size = length(u0)
    
    # Number of points in which to evaluate the solution or the BDF_weights
    r = 16
    N = ceil(Int, (tfinal-t0)/dt)
    Nr = ceil(Int, (N+1)/r)*r
    Q = ceil(Int, log2(Nr/r))-1
    NNr = 2^(Q+1)*r

    # Preallocation of some variables
    y = zeros(T, problem_size, N+1)
    fy = zeros(T, problem_size, N+1)
    zn = zeros(T, problem_size, NNr+1)

    # Evaluation of convolution and starting BDF_weights of the FLMM
    (omega, w, s) = BDF_weights(alpha, NNr+1)
    halpha = dt^alpha
    
    # Initializing solution and proces of computation
    mesh = collect(0:N)*dt
    y[:, 1] = u0
    temp = similar(u0)
    f(temp, u0, p, t0)
    fy[:, 1] = temp
    return BDFCache{iip, T}(prob, alg, mesh, u0, alpha, halpha, y, fy, zn,
                            p, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr,
                            omega, w, s, dt, reltol, abstol, maxiters, kwargs)
end

function SciMLBase.solve!(cache::BDFCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, halpha, y, fy, zn, p, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr, omega, w, s, dt, reltol, abstol, maxiters, kwargs = cache
    t0 = mesh[1]; tfinal = mesh[end]
    problem_size = length(u0)
    # generate jacobian of input function
    Jfdefun(t, u) = jacobian_of_fdefun(prob.f, t, u, p)

    BDF_first_approximations(cache, Jfdefun, t0)
    BDF_triangolo(cache, s+1, r-1, 0, N, Jfdefun, t0)
    
    # Main process of computation by means of the FFT algorithm
    nx0::Int = 0; ny0::Int = 0
    ff = zeros(T, 1, 2^(Q+2), 1)
    ff[1:2] = [0 2]
    for q = 0:Q
        L = 2^q
        BDF_disegna_blocchi(cache, L, ff, r, Nr, nx0+L*r, ny0, N, problem_size, Jfdefun, t0)
        ff[1:4*L] = [ff[1:2*L]; ff[1:2*L-1]; 4*L]
    end
    # Evaluation solution in TFINAL when TFINAL is not in the mesh
    if tfinal < mesh[N+1]
        c = (tfinal - mesh[N])/dt
        mesh[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end
    mesh = mesh[1:N+1]; y = y[:, 1:N+1]
    y = collect(Vector{eltype(y)}, eachcol(y))
    
    return DiffEqBase.build_solution(prob, alg, mesh, y)
end


function BDF_disegna_blocchi(cache, L, ff, r, Nr, nx0, ny0, N, problem_size, Jfdefun, t0)
    @unpack mesh, y, fy, zn, abstol, maxiters, s, w, omega, halpha, u0 = cache

    nxi::Int = copy(nx0); nxf::Int = copy(nx0 + L*r - 1)
    nyi::Int = copy(ny0); nyf::Int = copy(ny0 + L*r - 1)
    is::Int = 1
    T = eltype(cache)
    s_nxi = zeros(T, N)
    s_nxf = zeros(T, N)
    s_nyi = zeros(T, N)
    s_nyf = zeros(T, N)
    s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
    i_triangolo = 0;  stop = false
    while ~stop
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1)
        zn = BDF_quadrato(nxi, nxf, nyi, nyf, fy, zn, omega, problem_size)
        BDF_triangolo(cache, nxi, nxi+r-1, nxi, N, Jfdefun, t0)
        i_triangolo = i_triangolo + 1
        
        if ~stop
            if nxi+r-1 == nxf   # Il triangolo finisce dove finisce il quadrato -> si scende di livello
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is] + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            else # Il triangolo finisce prima del quadrato -> si fa un quadrato accanto
                nxi = nxi + r; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf = nyf + r
                is = is + 1
                s_nxi[is] = nxi; s_nxf[is] = nxf ; s_nyi[is] = nyi ; s_nyf[is] = nyf
            end
        end
        
    end
    return y, fy
end

function BDF_quadrato(nxi, nxf, nyi, nyf, fy, zn, omega, problem_size)
    coef_beg::Int = nxi-nyf; coef_end::Int = nxf-nyi+1
    funz_beg::Int = nyi+1; funz_end::Int = nyf+1
    vett_coef = omega[coef_beg+1:coef_end+1]
    vett_funz = [fy[:, funz_beg:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
    zzn = real(fast_conv(vett_coef, vett_funz))
    zn[:, nxi+1:nxf+1] = zn[:, nxi+1:nxf+1] + zzn[:, nxf-nyf:end-1]
    return zn
end

function BDF_triangolo(cache::BDFCache{iip, T}, nxi, nxf, j0, N, Jfdefun, t0) where{iip, T}
    @unpack prob, mesh, zn, abstol, maxiters, s, w, omega, halpha, u0, m_alpha, m_alpha_factorial, p = cache
    problem_size = length(u0)
    for n = nxi:min(N, nxf)
        n1::Int = n+1
        St = ABM_starting_term(mesh[n1], u0, m_alpha, t0, m_alpha_factorial)
        
        Phi = zeros(T, problem_size, 1)
        for j = 0:s
            Phi = Phi + w[j+1, n1]*cache.fy[:, j+1]
        end
        for j = j0:n-1
            Phi = Phi + omega[n-j+1]*cache.fy[:, j+1]
        end
        Phi_n = St + halpha*(zn[:, n1] + Phi)
        
        yn0 = cache.y[:, n]
        temp = zeros(T, problem_size)
        prob.f(temp, yn0, p, mesh[n1])
        fn0 = temp
        Jfn0 = Jf_vectorfield(mesh[n1], yn0, Jfdefun)
        Gn0 = yn0 - halpha*omega[1]*fn0 - Phi_n
        stop = false; it = 0
        yn1 = similar(yn0)
        fn1 = similar(yn0)
        while ~stop            
            JGn0 = zeros(T, problem_size, problem_size)+I - halpha*omega[1]*Jfn0
            yn1 = yn0 - JGn0\Gn0
            prob.f(fn1, yn1, p, mesh[n1])
            Gn1 = yn1 - halpha*omega[1]*fn1 - Phi_n
            it = it + 1
            
            stop = (norm(yn1-yn0, Inf) < abstol) || (norm(Gn1, Inf)<abstol)
            if it > maxiters && ~stop
                @warn "Non Convergence"
                stop = true
            end
            
            yn0 = yn1; Gn0 = Gn1
            if ~stop
                Jfn0 = Jf_vectorfield(mesh[n1], yn0, Jfdefun)
            end
            
        end
        cache.y[:, n1] = yn1
        cache.fy[:, n1] = fn1
    end
end

function BDF_first_approximations(cache::BDFCache{iip, T}, Jfdefun, t0) where {iip, T}
    @unpack prob, mesh, abstol, maxiters, s, halpha, omega, w, u0, m_alpha, m_alpha_factorial, p = cache
    m = length(u0)
    Im = zeros(m, m)+I; Ims = zeros(m*s, m*s)+I
    Y0 = zeros(T, s*m, 1); F0 = copy(Y0); B0 = copy(Y0)
    for j = 1 : s
        Y0[(j-1)*m+1:j*m, 1] = cache.y[:, 1]
        temp = zeros(length(cache.y[:, 1]))
        prob.f(temp, cache.y[:, 1], p, mesh[j+1])
        F0[(j-1)*m+1:j*m, 1] = temp
        St = ABM_starting_term(mesh[j+1], u0, m_alpha, t0, m_alpha_factorial)
        B0[(j-1)*m+1:j*m, 1] = St + halpha*(omega[j+1]+w[1, j+1])*cache.fy[:, 1]
    end
    W = zeros(T, s, s)
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
    JF = zeros(T, s*m, s*m)
    for j = 1:s
        JF[(j-1)*m+1:j*m, (j-1)*m+1:j*m] = Jf_vectorfield(mesh[j+1], cache.y[:, 1], Jfdefun)
    end
    stop = false; it = 0
    F1 = zeros(T, s*m, 1)
    while ~stop
        JG = Ims - W*JF
        global Y1 = Y0 - JG\G0
        
        for j = 1 : s
            temp = zeros(T, length(Y1[(j-1)*m+1:j*m, 1]))
            prob.f(temp, Y1[(j-1)*m+1:j*m, 1], p, mesh[j+1])
            F1[(j-1)*m+1:j*m, 1] = temp
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
                JF[(j-1)*m+1:j*m, (j-1)*m+1:j*m] = Jf_vectorfield(mesh[j+1], Y1[(j-1)*m+1:j*m, 1], Jfdefun)
            end
        end
        
    end
    for j = 1 : s
        cache.y[:, j+1] = Y1[(j-1)*m+1:j*m, 1]
        cache.fy[:, j+1] = F1[(j-1)*m+1:j*m, 1]
    end
end

function BDF_weights(alpha, N)
    # BDF-2 with generating function (2/3/(1-4x/3+x^2/3))^alpha
    omega = zeros(1, N+1)
    onethird = 1/3; fourthird = 4/3
    twothird_oneminusalpha = 2/3*(1-alpha)
    fourthird_oneminusalpha = 4/3*(1-alpha)
    omega[1] = 1; omega[2] = fourthird*alpha*omega[1]
    for n = 2:N
        omega[n+1] = (fourthird - fourthird_oneminusalpha/n)*omega[n] + (twothird_oneminusalpha/n - onethird)*omega[n-1]
    end
    omega = omega*((2/3)^(alpha))

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

Jf_vectorfield(t, y, Jfdefun) = Jfdefun(t, y)

function ABM_starting_term(t,y0, m_alpha, t0, m_alpha_factorial)
    ys = zeros(size(y0, 1), 1)
    for k = 1:Int64(m_alpha)
        ys = ys + (t-t0)^(k-1)/m_alpha_factorial[k]*y0[:, k]
    end
    return ys
end

function jacobian_of_fdefun(f, t, y, p)
    ForwardDiff.jacobian(y) do y # Thanks Bagge Carlson for helping me with this
    du = similar(y)
    f(du, y, p, t)
    du
    end
end