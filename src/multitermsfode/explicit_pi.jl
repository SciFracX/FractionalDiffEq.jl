@concrete mutable struct MTPIEXCache{T}
    prob
    alg
    mesh
    u0
    bet
    lam_rat_i
    gamma_val

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
    bn

    kwargs
end

Base.eltype(::MTPIEXCache{T}) where {T} = T

function SciMLBase.__init(
        prob::MultiTermsFODEProblem, alg::MTPIEX; dt = 0.0, abstol = 1e-6, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    (; parameters, orders, f, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    u0 = u0[:]'
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

    r::Int = 16
    N::Int = ceil(Int, (tfinal - t0) / dt)
    Nr::Int = ceil(Int, (N + 1) / r) * r
    Qr::Int = ceil(Int, log2((Nr) / r)) - 1
    NNr::Int = 2^(Qr + 1) * r

    y = zeros(Float64, problem_size, N + 1)
    fy = zeros(Float64, problem_size, N + 1)
    zn = zeros(Float64, problem_size, NNr + 1, orders_length)

    nvett = collect(0:(NNr + 1))
    bn = zeros(orders_length, NNr + 1)
    for i in 1:orders_length
        nbeta = nvett .^ bet[i]
        bn[i, :] = (nbeta[2:end] - nbeta[1:(end - 1)]) * dt^bet[i] / gamma(bet[i] + 1)
    end
    C = 0
    for i in 1:(orders_length - 1)
        C = C + lam_rat_i[i] * bn[i, 1]
    end

    mesh = collect(0:N) * dt
    y[:, 1] = u0[:, 1]
    fy[:, 1] .= f(u0[:], p, t0)

    return MTPIEXCache{T}(prob, alg, mesh, u0, bet, lam_rat_i, gamma_val,
        highest_order_parameter, highest_order_ceiling,
        other_orders_ceiling, y, fy, p, zn, r, N, Nr, Qr, NNr, bn, kwargs)
end

function SciMLBase.solve!(cache::MTPIEXCache{T}) where {T}
    (; prob, alg, mesh, y, r, N, Qr) = cache
    t0 = mesh[1]
    tfinal = mesh[end]
    MTPX_multiterms_triangolo(cache, 1, r - 1, t0)

    ff = zeros(1, 2^(Qr + 2))
    ff[1:2] = [0 2]
    card_ff = 2
    nx0 = 0
    nu0 = 0
    for qr in 0:Qr
        L = 2^qr
        MTPX_multiterms_disegna_blocchi(cache, L, ff, nx0 + L * r, nu0, t0)
        ff[1:(2 * card_ff)] = [ff[1:card_ff] ff[1:card_ff]]
        card_ff = 2 * card_ff
        ff[card_ff] = 4 * L
    end

    if tfinal < mesh[N + 1]
        c = (tfinal - mesh[N]) / dt
        mesh[N + 1] = tfinal
        y[:, N + 1] = (1 - c) * y[:, N] + c * y[:, N + 1]
    end

    mesh = mesh[1:(N + 1)]
    y = y[:, 1:(N + 1)]
    y = collect(Vector{eltype(y)}, eachcol(y))
    return DiffEqBase.build_solution(prob, alg, mesh, y)
end

function MTPX_multiterms_disegna_blocchi(
        cache::MTPIEXCache{T}, L, ff, nx0, nu0, t0) where {T}
    (; r, N, Nr) = cache
    nxi::Int = nx0
    nxf::Int = nx0 + L * r - 1
    nyi::Int = nu0
    nyf::Int = nu0 + L * r - 1
    is::Int = 1
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

        MTPX_multiterms_quadrato(cache, nxi, nxf, nyi, nyf)

        MTPX_multiterms_triangolo(cache, nxi, nxi + r - 1, t0)
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

function MTPX_multiterms_quadrato(cache::MTPIEXCache{T}, nxi, nxf, nyi, nyf) where {T}
    (; prob, bn) = cache
    orders_length = length(prob.orders)
    problem_size = size(prob.u0, 2)
    coef_beg = nxi - nyf
    coef_end = nxf - nyi + 1
    funz_beg = nyi + 1
    funz_end = nyf + 1

    for i in 1:orders_length
        vett_coef = bn[i, coef_beg:coef_end]
        if i < orders_length
            vett_funz = [cache.y[:, funz_beg:funz_end] zeros(problem_size,
                funz_end - funz_beg + 1)]
        else
            vett_funz = [cache.fy[:, funz_beg:funz_end] zeros(problem_size,
                funz_end - funz_beg + 1)]
        end
        zzn = real.(fast_conv(vett_coef, vett_funz))
        cache.zn[:, (nxi + 1):(nxf + 1), i] = cache.zn[:, (nxi + 1):(nxf + 1), i] +
                                              zzn[:, (nxf - nyf):(end - 1)]
    end
end

function MTPX_multiterms_triangolo(cache::MTPIEXCache{T}, nxi, nxf, t0) where {T}
    (; prob, mesh, u0, bet, lam_rat_i, gamma_val, highest_order_parameter, highest_order_ceiling, other_orders_ceiling, p, zn, N, bn) = cache
    problem_size = size(prob.u0, 2)
    orders_length = length(prob.orders)
    for n in nxi:min(N, nxf)
        Phi_n = MTPX_multiterms_starting_term(
            mesh[n + 1], t0, problem_size, u0, orders_length,
            highest_order_ceiling, other_orders_ceiling, bet, lam_rat_i, gamma_val)
        if nxi == 1
            j_beg = 0
        else
            j_beg = nxi
        end

        for i in 1:(orders_length - 1)
            temp = zn[:, n + 1, i]
            for j in j_beg:(n - 1)
                temp = temp + bn[i, n - j] * cache.y[:, j + 1]
            end
            Phi_n = Phi_n - lam_rat_i[i] * temp
        end
        temp = zn[:, n + 1, orders_length]
        for j in j_beg:(n - 1)
            temp = temp + bn[orders_length, n - j] * cache.fy[:, j + 1]
        end
        Phi_n = Phi_n + temp / highest_order_parameter

        cache.y[:, n + 1] = Phi_n
        cache.fy[:, n + 1] .= prob.f(cache.y[:, n + 1], p, mesh[n + 1])
    end
end

function MTPX_multiterms_starting_term(
        t, t0, problem_size, u0, orders_length, highest_order_ceiling,
        other_orders_ceiling, bet, lam_rat_i, gamma_val)
    ys = zeros(problem_size)

    for k in 0:(highest_order_ceiling - 1)
        ys = ys .+ (t - t0)^k ./ gamma_val[orders_length, k + 1] * u0[:, k + 1]
    end
    for i in 1:(orders_length - 1)
        temp = zeros(problem_size)
        for k in 0:(other_orders_ceiling[i] - 1)
            temp = temp .+ (t - t0)^(k + bet[i]) / gamma_val[i, k + 1] * u0[:, k + 1]
        end
        ys = ys + lam_rat_i[i] * temp
    end
    return ys'
end
