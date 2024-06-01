@concrete mutable struct MTPITrapCache{T}
    prob
    alg
    mesh
    u0
    bet
    lam_rat_i
    gamma_val

    J

    highest_order_parameter
    highest_order_ceiling
    other_orders_ceiling

    y
    fy
    p

    zn

    r
    N
    Nr
    Qr
    NNr
    C

    an
    a0
    abstol
    maxiters

    kwargs
end

Base.eltype(::MTPITrapCache{T}) where {T} = T

function SciMLBase.__init(prob::MultiTermsFODEProblem, alg::MTPITrap;
        dt = 0.0, abstol = 1e-6, maxiters = 100, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack parameters, orders, f, u0, tspan, p = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    u0 = u0[:]'

    # Generate the jacobian of the given function
    J_fun(x, y) = ForwardDiff.derivative(x -> f(y, p, x), x)

    orders_length = length(orders)
    orders = sort(orders)
    sorted_parameters_index = sortperm(orders, rev = true)
    parameters = parameters[sorted_parameters_index]
    highest_order = orders[end]
    other_orders = orders[1:(end - 1)]
    highest_order_parameter = parameters[end]
    lam_rat_i = parameters[1:(end - 1)] / highest_order_parameter
    highest_order_ceiling = ceil(Int64, highest_order)
    other_orders_ceiling = ceil.(Int64, orders[1:(end - 1)])
    bet = [highest_order .- other_orders; highest_order]

    gamma_val = zeros(orders_length, highest_order_ceiling)
    for i in 1:(orders_length - 1)
        k = collect(Int, 0:(other_orders_ceiling[i] - 1))
        gamma_val[i, k .+ 1] = gamma.(k .+ bet[i] .+ 1)
    end
    k = collect(0:(highest_order_ceiling - 1))
    gamma_val[orders_length, :] = factorial.(k)

    problem_size = size(u0, 1)

    r::Int = 64
    N::Int = ceil(Int, (tfinal - t0) / dt)
    Nr::Int = ceil(Int, (N + 1) / r) * r
    Qr::Int = ceil(Int, log2((Nr) / r)) - 1
    NNr::Int = 2^(Qr + 1) * r

    y = zeros(problem_size, N + 1)
    fy = zeros(problem_size, N + 1)
    zn = zeros(problem_size, NNr + 1, orders_length)

    nvett = collect(0:(NNr + 1))
    an = zeros(orders_length, NNr + 1)
    a0 = zeros(orders_length, NNr + 1)
    for i in 1:orders_length
        nbeta = nvett .^ bet[i]
        nbeta1 = nbeta .* nvett
        an[i, :] = [1; (nbeta1[1:(end - 2)] - 2 * nbeta1[2:(end - 1)] + nbeta1[3:end])] *
                   dt^bet[i] / gamma(bet[i] + 2)
        a0[i, :] = [0;
                    nbeta1[1:(end - 2)] -
                    nbeta[2:(end - 1)] .* (nvett[2:(end - 1)] .- bet[i] .- 1)] * dt^bet[i] /
                   gamma(bet[i] + 2)
    end
    C = 0
    for i in 1:(orders_length - 1)
        C = C + lam_rat_i[i] * an[i, 1]
    end

    mesh = collect(0:N) * dt
    y[:, 1] = u0[:, 1]
    fy[:, 1] .= f(u0[:, 1], p, t0)

    return MTPITrapCache{T}(prob, alg, mesh, u0, bet, lam_rat_i, gamma_val, J_fun,
        highest_order_parameter, highest_order_ceiling, other_orders_ceiling,
        y, fy, p, zn, r, N, Nr, Qr, NNr, C, an, a0, abstol, maxiters, kwargs)
end

function SciMLBase.solve!(cache::MTPITrapCache{T}) where {T}
    @unpack prob, alg, mesh, u0, bet, lam_rat_i, gamma_val, J, highest_order_parameter, highest_order_ceiling, other_orders_ceiling, y, fy, p, zn, r, N, Nr, Qr, NNr, C, an, a0, abstol, maxiters, kwargs = cache

    t0 = mesh[1]
    tfinal = mesh[end]
    MTPITrap_triangolo(cache, 1, r - 1, t0)

    ff = zeros(1, 2^(Qr + 2))
    ff[1:2] = [0 2]
    card_ff = 2
    nx0 = 0
    nu0 = 0
    for qr in 0:Qr
        L = 2^qr
        MTPITrap_disegna_blocchi(cache, L, ff, nx0 + L * r, nu0, t0)
        ff[1:(2 * card_ff)] = [ff[1:card_ff] ff[1:card_ff]]
        card_ff = 2 * card_ff
        ff[card_ff] = 4 * L
    end

    if tfinal < mesh[N + 1]
        c = (tfinal - mesh[N]) / dt
        mesh[N + 1] = tfinal
        cache.y[:, N + 1] = (1 - c) * cache.y[:, N] + c * cache.y[:, N + 1]
    end

    mesh = mesh[1:(N + 1)]
    u = cache.y[:, 1:(N + 1)]
    u = collect(Vector{eltype(u0)}, eachcol(u))

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end

function MTPITrap_disegna_blocchi(cache::MTPITrapCache{T}, L, ff, nx0, nu0, t0) where {T}
    @unpack prob, alg, mesh, u0, bet, lam_rat_i, gamma_val, J, highest_order_parameter, highest_order_ceiling, other_orders_ceiling, p, zn, r, N, Nr, Qr, NNr, C, an, a0, abstol, maxiters, kwargs = cache

    nxi::Int = nx0
    nxf::Int = nx0 + L * r - 1
    nyi::Int = nu0
    nyf::Int = nu0 + L * r - 1
    is = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[1] = nxi
    s_nxf[1] = nxf
    s_nyi[1] = nyi
    s_nyf[1] = nyf

    i_triangolo = 0
    stop = false
    while stop == false
        stop = (nxi + r - 1 == nx0 + L * r - 1) || (nxi + r - 1 >= Nr - 1)

        MTPITrap_quadrato(cache, nxi, nxf, nyi, nyf)

        MTPITrap_triangolo(cache, nxi, nxi + r - 1, t0)
        i_triangolo = i_triangolo + 1

        if stop == false
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

function MTPITrap_quadrato(cache::MTPITrapCache{T}, nxi, nxf, nyi, nyf) where {T}
    @unpack prob, u0, an = cache

    coef_beg = nxi - nyf
    coef_end = nxf - nyi + 1
    funz_beg = nyi + 1
    funz_end = nyf + 1

    orders_length = length(prob.orders)
    problem_size = size(u0, 1)

    for i in 1:orders_length
        vett_coef = an[i, coef_beg:coef_end]
        if nyi == 0
            if i < orders_length
                vett_funz = [zeros(problem_size, 1) cache.y[:, (funz_beg + 1):funz_end] zeros(problem_size,
                    funz_end - funz_beg + 1)]
            else
                vett_funz = [zeros(problem_size, 1) cache.fy[:, (funz_beg + 1):funz_end] zeros(problem_size,
                    funz_end - funz_beg + 1)]
            end
        else
            if i < orders_length
                vett_funz = [cache.y[:, funz_beg:funz_end] zeros(problem_size,
                    funz_end - funz_beg + 1)]
            else
                vett_funz = [cache.fy[:, funz_beg:funz_end] zeros(problem_size,
                    funz_end - funz_beg + 1)]
            end
        end
        zzn = real.(fast_conv(vett_coef, vett_funz))
        cache.zn[:, (nxi + 1):(nxf + 1), i] = cache.zn[:, (nxi + 1):(nxf + 1), i] +
                                              zzn[:, (nxf - nyf + 1):end]
    end
end

function MTPITrap_triangolo(cache::MTPITrapCache{T}, nxi, nxf, t0) where {T}
    @unpack prob, alg, mesh, u0, bet, lam_rat_i, gamma_val, J, highest_order_parameter, highest_order_ceiling, other_orders_ceiling, p, zn, r, N, Nr, Qr, NNr, C, an, a0, abstol, maxiters, kwargs = cache

    orders_length = length(prob.orders)
    problem_size = size(u0, 1)
    for n in nxi:min(N, nxf)
        St = MTPITrap_startingterm_multi(mesh[n + 1], t0, problem_size, u0, orders_length,
            highest_order_ceiling, other_orders_ceiling, bet, lam_rat_i, gamma_val)
        Phi_n = copy(St)

        for i in 1:(orders_length - 1)
            temp = a0[i, n + 1] * cache.y[:, 1] + zn[:, n + 1, i]
            for j in nxi:(n - 1)
                temp = temp + an[i, n - j + 1] * cache.y[:, j + 1]
            end
            Phi_n = Phi_n - lam_rat_i[i] * temp
        end
        temp = a0[orders_length, n + 1] * cache.fy[:, 1] + zn[:, n + 1, orders_length]
        for j in nxi:(n - 1)
            temp = temp + an[orders_length, n - j + 1] * cache.fy[:, j + 1]
        end
        Phi_n = Phi_n + temp / highest_order_parameter

        yn0 = cache.y[:, n]
        fn0 = prob.f(yn0, p, mesh[n + 1])
        Jfn0 = Jf_vectorfield(mesh[n + 1], yn0, J)
        Gn0 = (1 + C) * yn0 .- an[orders_length, 1] ./ highest_order_parameter * fn0 .-
              Phi_n
        stop = false
        it = 0

        while ~stop
            JGn0 = (1 + C) * (zeros(problem_size, problem_size) + I) .-
                   an[orders_length, 1] / highest_order_parameter * Jfn0
            global yn1 = yn0 - JGn0 \ Gn0
            global fn1 = prob.f(yn1, p, mesh[n + 1])
            Gn1 = (1 + C) * yn1 .- an[orders_length, 1] / highest_order_parameter * fn1 .-
                  Phi_n
            it = it + 1

            stop = (norm(yn1 - yn0, Inf) < abstol) || (norm(Gn1, Inf) < abstol)
            if it > maxiters && ~stop
                @warn "Non Convergence"
                stop = true
            end

            yn0 = yn1
            Gn0 = Gn1
            if ~stop
                Jfn0 = Jf_vectorfield(mesh[n + 1], yn0, J)
            end
        end

        cache.y[:, n + 1] = yn1
        cache.fy[:, n + 1] .= fn1
    end
end

function MTPITrap_startingterm_multi(
        mesh, t0, problem_size, u0, orders_length, highest_order_ceiling,
        other_orders_ceiling, bet, lam_rat_i, gamma_val)
    ys = zeros(problem_size)

    for k in 0:(highest_order_ceiling - 1)
        ys = ys .+ (mesh - t0)^k ./ gamma_val[orders_length, k + 1] * u0[:, k + 1]
    end
    for i in 1:(orders_length - 1)
        temp = zeros(problem_size)
        for k in 0:(other_orders_ceiling[i] - 1)
            temp = temp .+ (mesh - t0)^(k + bet[i]) / gamma_val[i, k + 1] * u0[:, k + 1]
        end
        ys = ys + lam_rat_i[i] * temp
    end
    return ys'
end
