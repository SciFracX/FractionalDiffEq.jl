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
    jac
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

Base.eltype(::BDFCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::BDF; dt = 0.0, reltol = 1e-6,
        abstol = 1e-6, maxiters = 1000, kwargs...)
    prob, iip = _is_need_convert!(prob)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
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

    # Number of points in which to evaluate the solution or the BDF_weights
    r = 16
    N = ceil(Int, (tfinal - t0) / dt)
    Nr = ceil(Int, (N + 1) / r) * r
    Q = ceil(Int, log2(Nr / r)) - 1
    NNr = 2^(Q + 1) * r

    # Preallocation of some variables
    y = [u0 for _ in 1:(N + 1)]
    fy = similar(y)
    zn = zeros(T, problem_size, NNr + 1)

    if prob.f.jac === nothing
        if iip
            jac = (df, u, p, t) -> begin
                _du = similar(u)
                prob.f(_du, u, p, t)
                _f = @closure (du, u) -> prob.f(du, u, p, t)
                ForwardDiff.jacobian!(df, _f, _du, u)
                return
            end
        else
            jac = (df, u, p, t) -> begin
                _du = prob.f(u, p, t)
                _f = @closure (du, u) -> (du .= prob.f(u, p, t))
                ForwardDiff.jacobian!(df, _f, _du, u)
                return
            end
        end
    else
        jac = prob.f.jac
    end

    # Evaluation of convolution and starting BDF_weights of the FLMM
    (omega, w, s) = BDF_weights(alpha, NNr + 1)
    halpha = dt^alpha

    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N) * dt
    y[1] .= high_order_prob ? u0[1, :] : u0
    temp = high_order_prob ? similar(u0[1, :]) : similar(u0)
    if iip
        prob.f(temp, u0, p, t0)
    else
        temp .= prob.f(u0, p, t0)
    end
    fy[1] = temp

    return BDFCache{iip, T}(prob, alg, mesh, u0, alpha, halpha, y, fy, zn, jac, prob.p,
        problem_size, m_alpha, m_alpha_factorial, r, N, Nr, Q, NNr, omega,
        w, s, dt, reltol, abstol, maxiters, high_order_prob, kwargs)
end

function SciMLBase.solve!(cache::BDFCache{iip, T}) where {iip, T}
    (; prob, alg, mesh, y, r, N, Q, s, dt) = cache

    BDF_first_approximations(cache)
    BDF_triangolo(cache, s + 1, r - 1, 0)

    # Main process of computation by means of the FFT algorithm
    nx0 = 0
    ny0 = 0
    ff = zeros(T, 1, 2^(Q + 2))
    ff[1:2] = [0; 2]
    for q in 0:Q
        L = 2^q
        BDF_disegna_blocchi(cache, L, ff, nx0 + L * r, ny0)
        ff[1:(4 * L)] = [ff[1:(2 * L)]; ff[1:(2 * L - 1)]; 4 * L]
    end

    return DiffEqBase.build_solution(prob, alg, mesh, cache.y)
end

function BDF_disegna_blocchi(
        cache::BDFCache{iip, T}, L::P, ff, nx0::P, ny0::P) where {P <: Integer, iip, T}
    (; r, Nr, N) = cache

    nxi::Int = nx0
    nxf::Int = nx0 + L * r - 1
    nyi::Int = ny0
    nyf::Int = ny0 + L * r - 1
    is = 1
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
        BDF_quadrato(cache, nxi, nxf, nyi, nyf)
        BDF_triangolo(cache, nxi, nxi + r - 1, nxi)
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

function BDF_quadrato(cache::BDFCache, nxi::P, nxf::P, nyi::P, nyf::P) where {P <: Integer}
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

function BDF_triangolo(
        cache::BDFCache{iip, T}, nxi::P, nxf::P, j0) where {P <: Integer, iip, T}
    (; prob, mesh, problem_size, zn, jac, N, abstol, maxiters, s, w, omega, halpha, p) = cache
    Jfn0 = Matrix{T}(undef, problem_size, problem_size)
    for n in nxi:min(N, nxf)
        n1 = n + 1
        St = ABM_starting_term(cache, mesh[n1])

        Phi = zeros(T, problem_size, 1)
        for j in 0:s
            Phi = Phi + w[j + 1, n1] * cache.fy[j + 1]
        end
        for j in j0:(n - 1)
            Phi = Phi + omega[n - j + 1] * cache.fy[j + 1]
        end
        Phi_n = St + halpha * (zn[:, n1] + Phi)
        yn0 = cache.y[n]
        fn0 = similar(yn0)
        prob.f(fn0, yn0, p, mesh[n1])
        @views jac(Jfn0, yn0, p, mesh[n1])
        Gn0 = yn0 - halpha * omega[1] * fn0 - Phi_n
        stop = false
        it = 0
        yn1 = similar(yn0)
        fn1 = similar(yn0)
        while !stop
            JGn0 = zeros(T, problem_size, problem_size) + I - halpha * omega[1] * Jfn0
            yn1 = yn0 - vec(JGn0 \ Gn0)
            prob.f(fn1, yn1, p, mesh[n1])
            Gn1 = yn1 - halpha * omega[1] * fn1 - Phi_n
            it = it + 1

            stop = (norm(yn1 - yn0, Inf) < abstol) || (norm(Gn1, Inf) < abstol)
            if it > maxiters && !stop
                @warn "Non Convergence"
                stop = true
            end

            yn0 = yn1
            Gn0 = Gn1
            if !stop
                @views jac(Jfn0, yn0, p, mesh[n1])
            end
        end
        cache.y[n1] = yn1
        cache.fy[n1] = fn1
    end
end

function BDF_first_approximations(cache::BDFCache{iip, T}) where {iip, T}
    (; prob, mesh, abstol, problem_size, maxiters, s, halpha, omega, w, jac, p) = cache

    Im = zeros(problem_size, problem_size) + I
    Ims = zeros(problem_size * s, problem_size * s) + I
    Y0 = VectorOfArray([cache.y[1] for _ in 1:s])
    F0 = similar(Y0)
    B0 = similar(Y0)
    for j in 1:s
        if iip
            prob.f(F0.u[j], cache.y[1], p, mesh[j + 1])
        else
            F0.u[j] .= prob.f(cache.y[1], p, mesh[j + 1])
        end
        St = ABM_starting_term(cache, mesh[j + 1])
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
    J_temp = Matrix{T}(undef, problem_size, problem_size)
    for j in 1:s
        if iip
            jac(J_temp, cache.y[1], p, mesh[j + 1])
        else
            J_temp .= jac(cache.y[1], p, mesh[j + 1])
        end
        JF[((j - 1) * problem_size + 1):(j * problem_size), ((j - 1) * problem_size + 1):(j * problem_size)] .= J_temp
    end
    stop = false
    it = 0
    F1 = similar(F0)
    Y1 = similar(Y0)
    while !stop
        JG = Ims - W * JF
        recursive_unflatten!(Y1, vec(Y0) - JG \ G0)

        for j in 1:s
            prob.f(F1.u[j], Y1.u[j], p, mesh[j + 1])
        end
        G1 = vec(Y1 - B0) - W * vec(F1)

        it = it + 1
        stop = (norm(Y1 - Y0, Inf) < abstol) || (norm(G1, Inf) < abstol)
        if it > maxiters && !stop
            @warn "Non Convergence"
            stop = 1
        end

        Y0 = Y1
        G0 = G1
        if ~stop
            for j in 1:s
                jac(J_temp, Y1.u[j], p, mesh[j+1])
                JF[((j - 1) * problem_size + 1):(j * problem_size), ((j - 1) * problem_size + 1):(j * problem_size)] .= J_temp
            end
        end
    end
    for j in 1:s
        cache.y[j + 1] = Y1.u[j]
        cache.fy[j + 1] = F1.u[j]
    end
end

function BDF_weights(alpha, N)
    # BDF-2 with generating function (2/3/(1-4x/3+x^2/3))^alpha
    omega = zeros(1, N + 1)
    onethird = 1 / 3
    fourthird = 4 / 3
    twothird_oneminusalpha = 2 / 3 * (1 - alpha)
    fourthird_oneminusalpha = 4 / 3 * (1 - alpha)
    omega[1] = 1
    omega[2] = fourthird * alpha * omega[1]
    for n in 2:N
        omega[n + 1] = (fourthird - fourthird_oneminusalpha / n) * omega[n] +
                       (twothird_oneminusalpha / n - onethird) * omega[n - 1]
    end
    omega = omega * ((2 / 3)^(alpha))

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
                nn_nu_alpha[i + 1, :] = zeros(N + 1)
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

Jf_vectorfield(t, y, Jfdefun) = Jfdefun(t, y)

function ABM_starting_term(cache::BDFCache{iip, T}, t) where {iip, T}
    (; u0, m_alpha, mesh, m_alpha_factorial, high_order_prob) = cache
    t0 = mesh[1]
    u0 = high_order_prob ? reshape(u0, 1, length(u0)) : u0
    ys = zeros(size(u0, 1))
    for k in 1:m_alpha
        ys = ys + (t - t0)^(k - 1) / m_alpha_factorial[k] * u0[:, k]
    end
    return ys
end

function jacobian_of_fdefun(f, t, y, p)
    ForwardDiff.jacobian(y) do y
        du = similar(y)
        f(du, y, p, t)
        du
    end
end

function _is_need_convert!(prob::FODEProblem)
    length(prob.u0) == 1 ? (_convert_single_term_to_vectorized_prob!(prob), true) : (prob, SciMLBase.isinplace(prob))
end

function _convert_single_term_to_vectorized_prob!(prob::FODEProblem)
    if SciMLBase.isinplace(prob)
        if isa(prob.u0, AbstractArray)
            new_prob = remake(prob; order = [prob.order])
        else
            new_prob = remake(prob; u0 = prob.u0, order = [prob.order])
        end
        return new_prob
    else
        function new_f(du, u, p, t)
            du[1] = prob.f(u[1], p, t)
        end
        new_fun = ODEFunction{true}(new_f) # make in-place
        new_prob = remake(prob; f = new_fun, u0 = [prob.u0], order = [prob.order])
        return new_prob
    end
end

@views function recursive_unflatten!(y::VectorOfArray, x::AbstractArray)
    i = 0
    for yᵢ in y
        copyto!(yᵢ, x[(i + 1):(i + length(yᵢ))])
        i += length(yᵢ)
    end
    return nothing
end