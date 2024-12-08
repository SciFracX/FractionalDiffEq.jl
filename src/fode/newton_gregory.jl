@concrete mutable struct NewtonGregoryCache{iip, T}
    prob
    alg
    mesh
    u0
    order
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

Base.eltype(::NewtonGregoryCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::NewtonGregory; dt = 0.0,
        reltol = 1e-6, abstol = 1e-6, maxiters = 1000, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    prob, iip = _is_need_convert!(prob)
    (; f, order, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    all(x -> x == order[1], order) ? nothing :
    throw(ArgumentError("BDF method is only for commensurate order FODE"))
    alpha = order[1] # commensurate ordre FODE

    m_alpha = ceil.(Int, alpha)
    m_alpha_factorial = factorial.(collect(0:(m_alpha - 1)))
    problem_size = length(order)
    u0_size = length(u0)
    high_order_prob = problem_size !== u0_size

    # Number of points in which to evaluate the solution or the weights
    r = 16
    N = ceil(Int, (tfinal - t0) / dt)
    Nr = ceil(Int, (N + 1) / r) * r
    Q = ceil(Int, log2((Nr) / r)) - 1
    NNr = 2^(Q + 1) * r

    # Preallocation of some variables
    y = [Vector{T}(undef, problem_size) for _ in 1:(N + 1)]
    fy = similar(y)
    zn = zeros(problem_size, NNr + 1)

    # generate jacobian of input function
    Jfdefun(t, u) = jacobian_of_fdefun(prob.f, t, u, p)

    # Evaluation of convolution and starting weights of the FLMM
    (omega, w, s) = NG_weights(alpha, NNr + 1)
    halpha = dt^alpha

    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N) * dt
    y[1] .= high_order_prob ? u0[1, :] : u0
    temp = high_order_prob ? similar(u0[1, :]) : similar(u0)
    f(temp, u0, p, t0)
    fy[1] = temp

    return NewtonGregoryCache{iip, T}(
        prob, alg, mesh, u0, alpha, halpha, y, fy, zn, Jfdefun, p,
        problem_size, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr, omega,
        w, s, dt, reltol, abstol, maxiters, high_order_prob, kwargs)
end
function SciMLBase.solve!(cache::NewtonGregoryCache{iip, T}) where {iip, T}
    (; prob, alg, mesh, u0, y, r, N, Q, s) = cache

    NG_first_approximations(cache)
    NG_triangolo(cache, s + 1, r - 1, 0)

    # Main process of computation by means of the FFT algorithm
    nx0::Int = 0
    ny0::Int = 0
    ff = zeros(1, 2^(Q + 2), 1)
    ff[1:2] = [0 2]
    for q in 0:Q
        L = Int64(2^q)
        NG_disegna_blocchi!(cache, L, ff, nx0 + L * r, ny0)
        ff[1:(4 * L)] = [ff[1:(2 * L)]; ff[1:(2 * L - 1)]; 4 * L]
    end

    return DiffEqBase.build_solution(prob, alg, mesh, y)
end

function NG_disegna_blocchi!(cache::NewtonGregoryCache{iip, T}, L::P, ff,
        nx0::P, ny0) where {P <: Integer, iip, T}
    (; r, Nr, N) = cache
    nxi::Int = nx0
    nxf::Int = nx0 + L * r - 1
    nyi::Int = ny0
    nyf::Int = ny0 + L * r - 1
    is::Int = 1
    s_nxi = Vector{T}(undef, N)
    s_nxf = similar(s_nxi)
    s_nyi = similar(s_nxi)
    s_nyf = similar(s_nxi)
    s_nxi[1] = nxi
    s_nxf[1] = nxf
    s_nyi[1] = nyi
    s_nyf[1] = nyf
    i_triangolo = 0
    stop = false
    while !stop
        stop = (nxi + r - 1 == nx0 + L * r - 1) || (nxi + r - 1 >= Nr - 1)
        NG_quadrato(cache, nxi, nxf, nyi, nyf)
        NG_triangolo(cache, nxi, nxi + r - 1, nxi)
        i_triangolo = i_triangolo + 1

        if !stop
            if nxi + r - 1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta * r
                nxi = s_nxf[is] + 1
                nxf = s_nxf[is] + Delta
                nyi = s_nxf[is] - Delta + 1
                nyf = s_nxf[is]
                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf
            else
                nxi = nxi + r
                nxf = nxi + r - 1
                nyi = nyf + 1
                nyf = nyf + r
                is = is + 1
                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf
            end
        end
    end
end

function NG_quadrato(cache::NewtonGregoryCache{iip, T}, nxi::P, nxf::P,
        nyi::P, nyf::P) where {P <: Integer, iip, T}
    (; problem_size, omega) = cache

    coef_beg = nxi - nyf
    coef_end = nxf - nyi + 1
    funz_beg = nyi + 1
    funz_end = nyf + 1
    vett_coef = omega[(coef_beg + 1):(coef_end + 1)]
    vett_funz = [reduce(hcat, cache.fy[funz_beg:funz_end]) zeros(problem_size, funz_end - funz_beg + 1)]
    zzn = real(fast_conv(vett_coef, vett_funz))
    cache.zn[:, (nxi + 1):(nxf + 1)] = cache.zn[:, (nxi + 1):(nxf + 1)] +
                                       zzn[:, (nxf - nyf):(end - 1)]
end

function NG_triangolo(
        cache::NewtonGregoryCache{iip, T}, nxi::P, nxf::P, j0) where {P <: Integer, iip, T}
    (; prob, mesh, problem_size, zn, Jfdefun, N, abstol, maxiters, s, w, omega, halpha, p) = cache
    for n in nxi:min(N, nxf)
        n1 = Int64(n + 1)
        St = NG_starting_term(cache, mesh[n1])

        Phi = zeros(problem_size, 1)
        for j in 0:s
            Phi = Phi + w[j + 1, n1] * cache.fy[j + 1]
        end
        for j in j0:(n - 1)
            Phi = Phi + omega[n - j + 1] * cache.fy[j + 1]
        end
        Phi_n = St + halpha * (zn[:, n1] + Phi)

        yn0 = cache.y[n]
        temp = zeros(length(yn0))
        prob.f(temp, yn0, p, mesh[n1])
        fn0 = temp
        Jfn0 = Jf_vectorfield(mesh[n1], yn0, Jfdefun)
        Gn0 = yn0 - halpha * omega[1] * fn0 - Phi_n
        stop = false
        it = 0
        yn1 = similar(yn0)
        fn1 = similar(yn0)
        while ~stop
            JGn0 = zeros(problem_size, problem_size) + I - halpha * omega[1] * Jfn0
            yn1 = yn0 - vec(JGn0 \ Gn0)
            prob.f(fn1, yn1, p, mesh[n1])
            Gn1 = yn1 - halpha * omega[1] * fn1 - Phi_n
            it = it + 1

            stop = (norm(yn1 - yn0, Inf) < abstol) || (norm(Gn1, Inf) < abstol)
            if it > maxiters && ~stop
                @warn "Non Convergence"
                stop = true
            end

            yn0 = yn1
            Gn0 = Gn1
            if ~stop
                Jfn0 = Jf_vectorfield(mesh[n1], yn0, Jfdefun)
            end
        end
        cache.y[n1] = yn1
        cache.fy[n1] = fn1
    end
end

function NG_first_approximations(cache::NewtonGregoryCache{iip, T}) where {iip, T}
    (; prob, mesh, abstol, problem_size, maxiters, s, halpha, omega, w, Jfdefun, p) = cache
    Im = zeros(problem_size, problem_size) + I
    Ims = zeros(problem_size * s, problem_size * s) + I
    Y0 = VectorOfArray([cache.y[1] for _ in 1:s])
    F0 = similar(Y0)
    B0 = similar(Y0)
    for j in 1:s
        prob.f(F0.u[j], cache.y[1], p, mesh[j + 1])
        St = NG_starting_term(cache, mesh[j + 1])
        B0.u[j] = St + halpha * (omega[j + 1] + w[1, j + 1]) * cache.fy[1]
    end
    W = zeros(T, s, s)
    for i in 1:s
        for j in 1:s
            if i >= j
                W[i, j] = omega[i - j + 1] + w[j + 1, i + 1]
            else
                W[i, j] = w[j + 1, i + 1]
            end
        end
    end
    W = halpha * kron(W, Im)
    G0 = vec(Y0 - B0) - W * vec(F0)
    JF = zeros(T, s * problem_size, s * problem_size)
    for j in 1:s
        JF[((j - 1) * problem_size + 1):(j * problem_size), ((j - 1) * problem_size + 1):(j * problem_size)] = Jf_vectorfield(
            mesh[j + 1], cache.y[1], Jfdefun)
    end
    stop = false
    it = 0
    F1 = similar(F0)
    Y1 = similar(Y0)
    while ~stop
        JG = Ims - W * JF
        recursive_unflatten!(Y1, vec(Y0) - JG \ G0)

        for j in 1:s
            prob.f(F1.u[j], Y1.u[j], p, mesh[j + 1])
        end
        G1 = vec(Y1 - B0) - W * vec(F1)

        it = it + 1

        stop = (norm(Y1 - Y0, Inf) < abstol) || (norm(G1, Inf) < abstol)
        if it > maxiters && ~stop
            @warn "Non Convergence"
            stop = 1
        end

        Y0 = Y1
        G0 = G1
        if ~stop
            for j in 1:s
                JF[((j - 1) * problem_size + 1):(j * problem_size), ((j - 1) * problem_size + 1):(j * problem_size)] = Jf_vectorfield(
                    mesh[j + 1], Y1.u[j], Jfdefun)
            end
        end
    end
    for j in 1:s
        cache.y[j + 1] = Y1.u[j]
        cache.fy[j + 1] = F1.u[j]
    end
end

function NG_weights(alpha, N)
    # Newton-Gregory formula with generating function (1-x)^(-alpha)*(1-alpha/2*(1-x))
    omega1 = zeros(1, N + 1)
    omega = copy(omega1)
    alphameno1 = alpha - 1
    omega1[1] = 1
    for n in 1:N
        omega1[n + 1] = (1 + alphameno1 / n) * omega1[n]
    end
    omega[1] = 1 - alpha / 2
    omega[2:(N + 1)] = (1 - alpha / 2) * omega1[2:(N + 1)] + alpha / 2 * omega1[1:N]

    k = floor(1 / abs(alpha))
    if abs(k - 1 / alpha) < 1.0e-12
        A = collect(0:k) * abs(alpha)
    else
        A = [collect(0:k) * abs(alpha); 1]
    end
    s = length(A) - 1
    # Generation of the matrix and the right hand--side vectors of the system
    nn = collect(0:N)
    V = zeros(s + 1, s + 1)
    jj_nu = zeros(s + 1, N + 1)
    nn_nu_alpha = copy(jj_nu)
    for i in 0:s
        nu = A[i + 1]
        V[i + 1, :] = collect(0:s) .^ nu
        jj_nu[i + 1, :] = nn .^ nu
        if alpha > 0
            nn_nu_alpha[i + 1, :] = gamma(nu + 1) / gamma(nu + 1 + alpha) *
                                    nn .^ (nu + alpha)
        else
            if i == 0
                nn_nu_alpha[i + 1, :] = zeros(1, N + 1)
            else
                nn_nu_alpha[i + 1, :] = gamma(nu + 1) / gamma(nu + 1 + alpha) *
                                        nn .^ (nu + alpha)
            end
        end
    end

    temp = fast_conv([omega zeros(size(omega))], [jj_nu zeros(size(jj_nu))])
    temp = real.(temp)
    b = nn_nu_alpha - temp[:, 1:(N + 1)]
    # Solution of the linear system with multiple right-hand side
    w = real.(V \ b)

    return omega, w, s
end

function NG_starting_term(cache::NewtonGregoryCache{iip, T}, t) where {iip, T}
    (; u0, m_alpha, mesh, m_alpha_factorial, high_order_prob) = cache
    t0 = mesh[1]
    u0 = high_order_prob ? reshape(u0, 1, length(u0)) : u0
    ys = zeros(size(u0, 1))
    for k in 1:Int64(m_alpha)
        ys = ys + (t - t0)^(k - 1) / m_alpha_factorial[k] * u0[:, k]
    end
    return ys
end
