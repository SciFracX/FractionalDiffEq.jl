@concrete mutable struct TrapezoidCache{iip, T}
    prob
    alg
    mesh
    u0
    order # temporary set as integer
    halpha
    y
    fy
    zn
    Jfdefun
    p
    problem_size

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
    high_order_prob
    
    kwargs
end

Base.eltype(::TrapezoidCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::Trapezoid;
                          dt = 0.0, reltol=1e-6, abstol=1e-6, maxiters=1000, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    prob = _is_need_convert!(prob)
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    all(x->x==order[1], order) ? nothing : throw(ArgumentError("BDF method is only for commensurate order FODE"))
    alpha = order[1] # commensurate ordre FODE
    (alpha > 1.0) && throw(ArgumentError("BDF method is only for order <= 1.0"))

    m_alpha = ceil.(Int, alpha)
    m_alpha_factorial = factorial.(collect(0:m_alpha-1))
    problem_size = length(order)
    u0_size = length(u0)
    high_order_prob = problem_size !== u0_size
    
    # Check compatibility size of the problem with size of the vector field
    
    # Number of points in which to evaluate the solution or the weights
    r = 16
    N = ceil(Int, (tfinal-t0)/dt)
    Nr = ceil(Int, (N+1)/r)*r
    Q = ceil(Int, log2((Nr)/r))-1
    NNr = 2^(Q+1)*r

    # Preallocation of some variables
    y = zeros(T, problem_size, N+1)
    fy = zeros(T, problem_size, N+1)
    zn = zeros(T, problem_size, NNr+1)

    # generate jacobian of input function
    Jfdefun(t, u) = jacobian_of_fdefun(prob.f, t, u, p)

    # Evaluation of convolution and starting weights of the FLMM
    (omega, w, s) = TrapWeights(alpha, NNr+1)
    halpha = dt^alpha
    
    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N)*dt
    y[:, 1] = high_order_prob ? u0[1, :] : u0
    temp = high_order_prob ? similar(u0[1, :]) : similar(u0)
    f(temp, u0, p, t0)
    fy[:, 1] = temp
    return TrapezoidCache{iip, T}(prob, alg, mesh, u0, alpha, halpha, y, fy, zn, Jfdefun,
                                  p, problem_size, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr,
                                  omega, w, s, dt, reltol, abstol, maxiters, high_order_prob, kwargs)
end

function SciMLBase.solve!(cache::TrapezoidCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, halpha, y, fy, zn, p, problem_size, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr, omega, w, s, dt, reltol, abstol, maxiters, kwargs = cache
    t0 = mesh[1]; tfinal = mesh[end]

    TrapFirstApproximations(cache)
    TrapTriangolo(cache, s+1, r-1, 0)
    
    # Main process of computation by means of the FFT algorithm
    nx0 = 0; ny0 = 0
    ff = zeros(T, 1, 2^(Q+2), 1)
    ff[1:2] = [0 2]
    for q = 0:Q
        L = 2^q
        TrapDisegnaBlocchi(cache, L, ff, nx0+L*r, ny0)
        ff[1:4*L] = [ff[1:2*L]; ff[1:2*L-1]; 4*L]
    end
    # Evaluation solution in TFINAL when TFINAL is not in the mesh
    if tfinal < mesh[N+1]
        c = (tfinal - mesh[N])/dt
        mesh[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end
    mesh = mesh[1:N+1]; y = y[:, 1:N+1]
    y = collect(Vector{eltype(u0)}, eachcol(y))

    return DiffEqBase.build_solution(prob, alg, mesh, y)
end


function TrapDisegnaBlocchi(cache::TrapezoidCache{iip, T}, L::P, ff, nx0::P, ny0::P) where {P <: Integer, iip, T}
    @unpack mesh, y, fy, zn, abstol, maxiters, r, Nr, N, Jfdefun, s, w, omega, halpha, u0 = cache

    nxi::Int = copy(nx0); nxf::Int = copy(nx0 + L*r - 1)
    nyi::Int = copy(ny0); nyf::Int = copy(ny0 + L*r - 1)
    is = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[is] = nxi; s_nxf[is] = nxf ; s_nyi[is] = nyi; s_nyf[is] = nyf
    i_triangolo = 0;  stop = false
    while ~stop
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1)
        TrapQuadrato(cache, nxi, nxf, nyi, nyf)
        TrapTriangolo(cache, nxi, nxi+r-1, nxi)
        i_triangolo = i_triangolo + 1
        
        if ~stop
            if nxi+r-1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is] + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            else
                nxi = nxi + r; nxf = nxi+r-1; nyi = nyf+1; nyf = nyf+r
                is = is + 1
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            end
        end
        
    end
end

function TrapQuadrato(cache::TrapezoidCache{iip, T}, nxi::P, nxf::P, nyi::P, nyf::P) where {P <: Integer, iip, T}
    @unpack problem_size, omega = cache

    coef_beg = nxi-nyf; coef_end = nxf-nyi+1
    funz_beg = nyi+1; funz_end = nyf+1
    vett_coef = omega[coef_beg+1:coef_end+1]
    vett_funz = [cache.fy[:, funz_beg:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
    zzn = real(fast_conv(vett_coef, vett_funz))
    cache.zn[:, nxi+1:nxf+1] = cache.zn[:, nxi+1:nxf+1] + zzn[:, nxf-nyf:end-1]
end

function TrapTriangolo(cache::TrapezoidCache{iip, T}, nxi::P, nxf::P, j0) where {P <: Integer, iip, T}
    @unpack prob, mesh, problem_size, zn, Jfdefun, N, abstol, maxiters, s, w, omega, halpha, u0, m_alpha, m_alpha_factorial, p = cache
    t0 = mesh[1]
    for n = nxi:min(N, nxf)
        n1::Int = n+1
        St = TrapStartingTerm(cache, mesh[n1])
        
        Phi = zeros(problem_size, 1)
        for j = 0:s
            Phi = Phi + w[j+1, n1]*cache.fy[:, j+1]
        end
        for j = j0:n-1
            Phi = Phi + omega[n-j+1]*cache.fy[:, j+1]
        end
        Phi_n = St + halpha*(zn[:, n1] + Phi)
        
        yn0 = cache.y[:, n]
        temp = zeros(length(yn0))
        prob.f(temp, yn0, p, mesh[n1])
        fn0 = temp
        Jfn0 = Jf_vectorfield(mesh[n1], yn0, Jfdefun)
        Gn0 = yn0 - halpha*omega[1]*fn0 - Phi_n
        stop = false; it = 0
        yn1 = similar(yn0)
        fn1 = similar(yn0)
        while ~stop
            JGn0 = zeros(problem_size, problem_size)+I - halpha*omega[1]*Jfn0
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

function TrapFirstApproximations(cache::TrapezoidCache{iip, T}) where {iip, T}
    @unpack prob, mesh, abstol, problem_size, maxiters, s, halpha, omega, w, Jfdefun, u0, m_alpha, m_alpha_factorial, p = cache
    Im = zeros(problem_size, problem_size)+I ; Ims = zeros(problem_size*s, problem_size*s)+I
    Y0 = zeros(s*problem_size, 1); F0 = copy(Y0); B0 = copy(Y0)
    for j = 1 : s
        Y0[(j-1)*problem_size+1:j*problem_size, 1] = cache.y[:, 1]
        temp = zeros(length(cache.y[:, 1]))
        prob.f(temp, cache.y[:, 1], p, mesh[j+1])
        F0[(j-1)*problem_size+1:j*problem_size, 1] = temp
        St = TrapStartingTerm(cache, mesh[j+1])
        B0[(j-1)*problem_size+1:j*problem_size, 1] = St + halpha*(omega[j+1]+w[1, j+1])*cache.fy[:, 1]
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
    JF = zeros(s*problem_size, s*problem_size)
    for j = 1:s
        JF[(j-1)*problem_size+1:j*problem_size, (j-1)*problem_size+1:j*problem_size] = Jf_vectorfield(mesh[j+1], cache.y[:, 1], Jfdefun)
    end
    stop = false; it = 0
    F1 = zeros(s*problem_size, 1)
    while ~stop
        JG = Ims - W*JF
        global Y1 = Y0 - JG\G0
        
        for j in 1:s
            temp = zeros(length(Y1[(j-1)*problem_size+1:j*problem_size, 1]))
            prob.f(temp, Y1[(j-1)*problem_size+1:j*problem_size, 1], p, mesh[j+1])
            F1[(j-1)*problem_size+1:j*problem_size, 1] = temp
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
                JF[(j-1)*problem_size+1:j*problem_size, (j-1)*problem_size+1:j*problem_size] = Jf_vectorfield(mesh[j+1], Y1[(j-1)*problem_size+1:j*problem_size, 1], Jfdefun)
            end
        end
    end
    for j = 1 : s
        cache.y[:, j+1] = Y1[(j-1)*problem_size+1:j*problem_size, 1]
        cache.fy[:, j+1] = F1[(j-1)*problem_size+1:j*problem_size, 1]
    end
end

function TrapWeights(alpha, N)
    # Trapezoid method with generating function ((1+x)/2/(1-x))^alpha
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

function TrapStartingTerm(cache::TrapezoidCache{iip, T}, t) where {iip, T}
    @unpack u0, m_alpha, mesh, m_alpha_factorial, high_order_prob = cache
    t0 = mesh[1]
    u0 = high_order_prob ? reshape(u0, 1, length(u0)) : u0
    ys = zeros(size(u0, 1), 1)
    for k = 1:m_alpha
        ys = ys + (t-t0)^(k-1)/m_alpha_factorial[k]*u0[:, k]
    end
    return ys
end