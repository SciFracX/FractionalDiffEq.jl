@concrete mutable struct PIEXCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    m_alpha
    m_alpha_factorial
    y
    fy
    p
    problem_size
    zn

    r
    N
    Nr
    Qr
    NNr

    bn

    halpha1
    mu
    abstol
    index_fft
    bn_fft
    high_order_prob

    kwargs
end

Base.eltype(::PIEXCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::PIEX; dt = 0.0, abstol = 1e-6, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    prob, iip = _is_need_convert!(prob)
    (; f, order, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    mu = 1

    alpha_length = length(order)
    order = (alpha_length == 1) ? order : order[:]
    problem_size = length(order)
    u0_size = length(u0)
    high_order_prob = problem_size !== u0_size

    m_alpha = ceil.(Int, order)
    m_alpha_factorial = zeros(alpha_length, maximum(m_alpha))
    for i in 1:alpha_length
        for j in 0:(m_alpha[i] - 1)
            m_alpha_factorial[i, j + 1] = factorial(j)
        end
    end

    r = 16
    N = ceil(Int64, (tfinal - t0) / dt)
    Nr = ceil(Int64, (N + 1) / r) * r
    Qr = ceil(Int64, log2(Nr / r)) - 1
    NNr = 2^(Qr + 1) * r

    # Preallocation of some variables
    y = zeros(T, problem_size, N + 1)
    fy = zeros(T, problem_size, N + 1)
    zn = zeros(problem_size, NNr + 1)

    # Evaluation of coefficients of the PECE method
    nvett = 0:(NNr + 1) |> collect
    bn = zeros(alpha_length, NNr + 1)
    for i_alpha in 1:alpha_length
        find_alpha = Float64[]
        if alpha_length == 1
            nothing
        else
            if order[i_alpha] == order[1:(i_alpha - 1)]
                push!(find_alpha, i_alpha)
            end
        end

        if isempty(find_alpha) == false
            bn[i_alpha, :] = bn[find_alpha[1], :]
        elseif abs(order[i_alpha] - 1) < 1e-14
            bn[i_alpha, :] = [1; zeros(NNr)]
        else
            nalpha = nvett .^ order[i_alpha]
            bn[i_alpha, :] = nalpha[2:end] - nalpha[1:(end - 1)]
        end
    end
    halpha1 = dt .^ order ./ gamma.(order .+ 1)

    if Qr >= 0
        index_fft = zeros(Int, 2, Qr + 1)
        for l in 1:(Qr + 1)
            if l == 1
                index_fft[1, l] = 1
                index_fft[2, l] = r * 2
            else
                index_fft[1, l] = index_fft[2, l - 1] + 1
                index_fft[2, l] = index_fft[2, l - 1] + 2^l * r
            end
        end

        bn_fft = zeros(Complex, alpha_length, index_fft[2, end])
        for l in 1:(Qr + 1)
            coef_end = 2^l * r
            for i_alpha in 1:alpha_length
                find_alpha = Float64[]
                if alpha_length == 1
                    nothing
                else
                    if order[i_alpha] == order[1:(i_alpha - 1)]
                        push!(find_alpha, i_alpha)
                    end
                end
                if isempty(find_alpha) == false
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = bn_fft[
                        find_alpha[1], index_fft[1, l]:index_fft[2, l]]
                else
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(
                        bn[i_alpha, 1:coef_end], coef_end)
                end
            end
        end
    else
        index_fft = 0
        bn_fft = 0
    end

    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N) * dt
    y[:, 1] = high_order_prob ? u0[1, :] : u0
    temp = high_order_prob ? similar(u0[1, :]) : similar(u0)
    f(temp, u0, p, t0)
    fy[:, 1] = temp
    return PIEXCache{iip, T}(prob, alg, mesh, u0, order, m_alpha, m_alpha_factorial, y,
        fy, p, problem_size, zn, r, N, Nr, Qr, NNr, bn, halpha1,
        mu, abstol, index_fft, bn_fft, high_order_prob, kwargs)
end
function SciMLBase.solve!(cache::PIEXCache{iip, T}) where {iip, T}
    (; prob, alg, mesh, u0, y, r, N, Qr) = cache
    tfinal = mesh[end]
    PIEX_triangolo(cache, 1, r - 1)

    # Main process of computation by means of the FFT algorithm
    ff = zeros(1, 2^(Qr + 2))
    ff[1:2] = [0; 2]
    card_ff = 2
    nx0 = 0
    ny0 = 0
    for qr in 0:Qr
        L = 2^qr
        PIEX_disegna_blocchi(cache, L, ff, nx0 + L * r, ny0)
        ff[1:(2 * card_ff)] = [ff[1:card_ff]; ff[1:card_ff]]
        card_ff = 2 * card_ff
        ff[card_ff] = 4 * L
    end

    # Evaluation solution in tfinal when tfinal is not in the mesh
    if tfinal < mesh[N + 1]
        c = (tfinal - mesh[N]) / dt
        mesh[N + 1] = tfinal
        y[:, N + 1] = (1 - c) * y[:, N] + c * y[:, N + 1]
    end

    mesh = mesh[1:(N + 1)]
    y = cache.y[:, 1:(N + 1)]
    u = collect(Vector{eltype(u0)}, eachcol(y))

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end

function PIEX_disegna_blocchi(
        cache::PIEXCache{iip, T}, L::P, ff, nx0::P, ny0::P) where {P <: Integer, iip, T}
    (; N, r, Nr) = cache

    nxi::Int = nx0
    nxf::Int = nx0 + L * r - 1
    nyi::Int = ny0
    nyf::Int = ny0 + L * r - 1
    is = 1
    s_nxi = zeros(T, N)
    s_nxf = zeros(T, N)
    s_nyi = zeros(T, N)
    s_nyf = zeros(T, N)
    s_nxi[is] = nxi
    s_nxf[is] = nxf
    s_nyi[is] = nyi
    s_nyf[is] = nyf

    i_triangolo = 0
    stop = false
    while stop == false
        stop = (nxi + r - 1 == nx0 + L * r - 1) || (nxi + r - 1 >= Nr - 1)

        PIEX_quadrato(cache, nxi, nxf, nyi, nyf)
        PIEX_triangolo(cache, nxi, nxi + r - 1)
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

function PIEX_quadrato(cache::PIEXCache{iip, T}, nxi::P, nxf::P,
        nyi::P, nyf::P) where {P <: Integer, iip, T}
    (; prob, r, N, problem_size, index_fft, bn_fft) = cache
    coef_end = nxf - nyi + 1
    alpha_length = length(prob.order)
    i_fft::Int = log2(coef_end / r)
    funz_beg = nyi + 1
    funz_end = nyf + 1
    Nnxf = min(N, nxf)

    # Evaluation convolution segment for the predictor
    vett_funz = cache.fy[:, funz_beg:funz_end]
    vett_funz_fft = rowfft(vett_funz, coef_end)
    zzn = zeros(problem_size, coef_end)
    for i in 1:problem_size
        i_alpha::Int = min(alpha_length, i)
        if abs(prob.order[i_alpha] - 1) > 1e-14
            Z = bn_fft[i_alpha, index_fft[1, i_fft]:index_fft[2, i_fft]] .*
                vett_funz_fft[i, :]
            zzn[i, :] = real.(ourifft(Z, coef_end))
        end
    end
    zzn = zzn[:, (nxf - nyf):(end - 1)]
    cache.zn[:, (nxi + 1):(Nnxf + 1)] = cache.zn[:, (nxi + 1):(Nnxf + 1)] +
                                        zzn[:, 1:(Nnxf - nxi + 1)]
end

function PIEX_triangolo(
        cache::PIEXCache{iip, T}, nxi::P, nxf::P) where {P <: Integer, iip, T}
    (; prob, mesh, order, p, problem_size, N, bn, halpha1) = cache

    alpha_length = length(order)
    for n in nxi:min(N, nxf)
        St = PIEX_system_starting_term(cache, mesh[n + 1])
        # Evaluation of the predictor
        Phi = zeros(problem_size, 1)
        if nxi == 1 # Case of the first triangle
            j_beg = 0
        else # Case of any triangle but not the first
            j_beg = nxi
        end
        for j in j_beg:(n - 1)
            Phi = Phi + bn[1:alpha_length, n - j] .* cache.fy[:, j + 1]
        end

        i_alpha_1 = findall(alpha -> abs(alpha - 1) < 1e-14, order)
        Phi[i_alpha_1] = cache.fy[i_alpha_1, n]
        St[i_alpha_1] = cache.y[i_alpha_1, n]

        cache.y[:, n + 1] = St + halpha1 .* (cache.zn[:, n + 1] + Phi)
        temp = zeros(length(cache.y[:, n + 1]))
        prob.f(temp, cache.y[:, n + 1], p, mesh[n + 1])
        cache.fy[:, n + 1] = temp
    end
end

function PIEX_system_starting_term(cache::PIEXCache{iip, T}, t) where {iip, T}
    (; mesh, u0, m_alpha, m_alpha_factorial, high_order_prob) = cache
    t0 = mesh[1]
    u0 = high_order_prob ? reshape(u0, 1, length(u0)) : u0
    ys = zeros(size(u0, 1), 1)
    for k in 1:maximum(m_alpha)
        if length(m_alpha) == 1
            ys = ys .+ (t - t0)^(k - 1) / m_alpha_factorial[k] * u0[:, k]
        else
            i_alpha = findall(x -> x >= k, m_alpha)
            ys[i_alpha] = ys[i_alpha] +
                          (t - t0)^(k - 1) * u0[i_alpha, k] ./ m_alpha_factorial[i_alpha, k]
        end
    end
    return ys
end
