@concrete mutable struct PECECache{iip, T}
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

    zn_pred
    zn_corr

    r
    N
    Nr
    Qr
    NNr

    an
    bn
    a0
    halpha1
    halpha2
    mu
    abstol
    index_fft
    an_fft
    bn_fft

    kwargs
end

Base.eltype(::PECECache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::PECE; dt = 0.0, abstol = 1e-6, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    mu = 1
    order = order[:]

    # Check compatibility size of the problem with number of fractional orders
    alpha_length = length(order)
    problem_size = length(u0)

    m_alpha = ceil.(Int, order)
    m_alpha_factorial = zeros(alpha_length, maximum(m_alpha))
    for i = 1 : alpha_length
        for j = 0 : m_alpha[i]-1
            m_alpha_factorial[i, j+1] = factorial(j)
        end
    end

    f_temp = zeros(size(u0[:, 1]))
    f(f_temp, u0, p, t0)

    r::Int = 16
    N::Int = ceil(Int64, (tfinal-t0)/dt)
    Nr::Int = ceil(Int64, (N+1)/r)*r
    Qr::Int = ceil(Int64, log2(Nr/r)) - 1
    NNr::Int = 2^(Qr+1)*r

    # Preallocation of some variables
    y = zeros(problem_size, N+1)
    fy = zeros(problem_size, N+1)
    zn_pred = zeros(problem_size, NNr+1)
    mu > 0 ?  (zn_corr = zeros(problem_size, NNr+1)) : (zn_corr = 0)

    # Evaluation of coefficients of the PECE method
    nvett = 0:NNr+1 |> collect
    bn = zeros(alpha_length, NNr+1); an = copy(bn); a0 = copy(bn)
    for i_alpha = 1:alpha_length
        find_alpha = Float64[]
        if order[i_alpha] == order[1:i_alpha-1]
            push!(find_alpha, i_alpha)
        end

        if isempty(find_alpha) == false
            bn[i_alpha, :] = bn[find_alpha[1], :]
            an[i_alpha, :] = an[find_alpha[1], :]
            a0[i_alpha, :] = a0[find_alpha[1], :]
        else
            nalpha = nvett.^order[i_alpha]
            nalpha1 = nalpha.*nvett
            bn[i_alpha, :] = nalpha[2:end] - nalpha[1:end-1]
            an[i_alpha, :] = [1; (nalpha1[1:end-2] - 2*nalpha1[2:end-1] + nalpha1[3:end]) ]
            a0[i_alpha, :] = [0; nalpha1[1:end-2]-nalpha[2:end-1].*(nvett[2:end-1].-order[i_alpha].-1)]
        end
    end
    halpha1 = dt.^order./gamma.(order.+1)
    halpha2 = dt.^order./gamma.(order.+2)

    # Evaluation of FFT of coefficients of the PECE method
    index_fft = zeros(Int, 2, Qr+1)
    if Qr >= 0
        for l = 1 : Qr+1
            if l == 1
                index_fft[1, l] = 1; index_fft[2, l] = r*2
            else
                index_fft[1, l] = index_fft[2, l-1]+1; index_fft[2, l] = index_fft[2, l-1]+2^l*r
            end
        end
        
        bn_fft = zeros(Complex, alpha_length, index_fft[2, Qr+1]); an_fft = copy(bn_fft)
        for l = 1 : Qr+1
            coef_end = 2^l*r
            for i_alpha = 1 : alpha_length
                find_alpha = Float64[]
                if order[i_alpha] == order[1:i_alpha-1]
                    push!(find_alpha, i_alpha)
                end
                if isempty(find_alpha) == false
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = bn_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = an_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
                else
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(bn[i_alpha, 1:coef_end], coef_end)
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(an[i_alpha, 1:coef_end], coef_end)
                end
            end
        end
    else
        bn_fft = 0; an_fft = 0
    end

    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N)*dt
    y[:, 1] = u0[:, 1]
    fy[:, 1] = f_temp
    return PECECache{iip, T}(prob, alg, mesh, u0, order, m_alpha, m_alpha_factorial, y, fy, p,
                             zn_pred, zn_corr,
                             r, N, Nr, Qr, NNr,
                             an, bn, a0, halpha1, halpha2,
                             mu, abstol, index_fft, an_fft, bn_fft, kwargs)
end
function SciMLBase.solve!(cache::PECECache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, y, fy, zn_pred, zn_corr, r, N, Nr, Qr, NNr, an, bn, a0, halpha1, halpha2, mu, abstol, index_fft, an_fft, bn_fft, kwargs = cache
    t0 = mesh[1]; tfinal = mesh[N+1]
    ABM_triangolo(cache, 1, r-1, t0)

    # Main process of computation by means of the FFT algorithm
    ff = zeros(T, 1, 2^(Qr+2)); ff[1:2] = [0; 2] ; card_ff = 2
    nx0::Int = 0; ny0::Int = 0
    for qr = 0 : Qr
        L = 2^qr 
        DisegnaBlocchi(cache, L, ff, r, Nr, nx0+L*r, ny0, t0)
        ff[1:2*card_ff] = [ff[1:card_ff]; ff[1:card_ff]] 
        card_ff = 2*card_ff
        ff[card_ff] = 4*L
    end

    # Evaluation solution in T when T is not in the mesh
    if tfinal < mesh[N+1]
        c = (tfinal - mesh[N])/dt
        mesh[N+1] = tfinal
        cache.y[:, N+1] = (1-c)*cache.y[:, N] + c*cache.y[:, N+1]
    end

    mesh = mesh[1:N+1]; u = cache.y[:, 1:N+1]
    u = collect(Vector{eltype(u0)}, eachcol(u))

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end


function DisegnaBlocchi(cache::PECECache{iip, T}, L, ff, r, Nr, nx0, ny0, t0) where {iip, T}
    @unpack N = cache
    nxi::Int = nx0; nxf::Int = nx0 + L*r - 1
    nyi::Int = ny0; nyf::Int = ny0 + L*r - 1
    is::Int = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf

    i_triangolo = 0; stop = false
    while stop == false
        
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1)
        
        ABM_quadrato(cache, nxi, nxf, nyi, nyf)

        ABM_triangolo(cache, nxi, nxi+r-1, t0)
        i_triangolo = i_triangolo + 1

        if stop == false
            if nxi+r-1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is]  + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf ;
            else # Il triangolo finisce prima del quadrato -> si fa un quadrato accanto
                nxi = nxi + r ; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf = nyf + r  ;
                is = is + 1 ;
                s_nxi[is] = nxi ; s_nxf[is] = nxf ; s_nyi[is] = nyi ; s_nyf[is] = nyf ;
            end
        end
        
    end
end

function ABM_quadrato(cache::PECECache{iip, T}, nxi, nxf, nyi, nyf) where {iip, T}
    @unpack prob, mesh, r, N, Nr, Qr, NNr, an, bn, a0, halpha1, halpha2, mu, abstol, index_fft, an_fft, bn_fft = cache
    problem_size = length(prob.u0)
    alpha_length = length(prob.order)
    coef_end::Int = nxf-nyi+1
    i_fft::Int = log2(coef_end/r)
    funz_beg::Int = nyi+1
    funz_end::Int = nyf+1
    Nnxf::Int = min(N, nxf)

    # Evaluation convolution segment for the predictor
    vett_funz = cache.fy[:, funz_beg:funz_end]
    vett_funz_fft = rowfft(vett_funz, coef_end)
    zzn_pred = zeros(problem_size, coef_end)
    for i = 1:problem_size
        i_alpha = min(alpha_length,i)
        Z = bn_fft[i_alpha, index_fft[1, i_fft]:index_fft[2, i_fft]].*vett_funz_fft[i, :]
        zzn_pred[i, :] = real.(ourifft(Z, coef_end))
    end
    zzn_pred = zzn_pred[:, nxf-nyf:end-1]
    cache.zn_pred[:, nxi+1:Nnxf+1] = cache.zn_pred[:, nxi+1:Nnxf+1] + zzn_pred[:, 1:Nnxf-nxi+1]

    # Evaluation convolution segment for the corrector
    if mu > 0
        if nyi == 0 # Evaluation of the lowest square
            vett_funz = [zeros(problem_size, 1) cache.fy[:, funz_beg+1:funz_end]]
            vett_funz_fft = rowfft(vett_funz, coef_end)
        end
        zzn_corr = zeros(problem_size, coef_end)
        for i=1:problem_size
            i_alpha = min(alpha_length,i)
            Z = an_fft[i_alpha, index_fft[1, i_fft]:index_fft[2, i_fft]].*vett_funz_fft[i, :]
            zzn_corr[i, :] = real.(ourifft(Z, coef_end))
        end
        zzn_corr = zzn_corr[:, nxf-nyf+1:end]
        cache.zn_corr[:, nxi+1:Nnxf+1] = cache.zn_corr[:, nxi+1:Nnxf+1] + zzn_corr[:, 1:Nnxf-nxi+1]
    else
        cache.zn_corr = 0
    end
end



function ABM_triangolo(cache::PECECache{iip, T}, nxi, nxf, t0) where {iip, T}
    @unpack prob, mesh, u0, order, m_alpha, m_alpha_factorial, p, zn_pred, zn_corr, N, an, bn, a0, halpha1, halpha2, mu, abstol, index_fft, an_fft, bn_fft = cache
    alpha_length = length(order)
    problem_size = length(u0)

    for n = nxi:min(N, nxf)
        # Evaluation of the predictor
        Phi = zeros(T, problem_size, 1)
        if nxi == 1 # Case of the first triangle
            j_beg::Int = 0
        else # Case of any triangle but not the first
            j_beg = nxi
        end
        for j = j_beg:n-1
            Phi = Phi + bn[1:alpha_length,n-j].*cache.fy[:, j+1]
        end
        St = starting_term(mesh[n+1], u0, m_alpha, t0, m_alpha_factorial)
        y_pred = St + halpha1.*(zn_pred[:, n+1] + Phi)
        f_pred = zeros(length(y_pred))
        prob.f(f_pred, y_pred, p, mesh[n+1])
        
        # Evaluation of the corrector
        if mu == 0
            cache.y[:, n+1] = y_pred
            cache.fy[:, n+1] = f_pred
        else
            j_beg = nxi
            Phi = zeros(problem_size, 1)
            for j = j_beg : n-1
                Phi = Phi + an[1:alpha_length, n-j+1].*cache.fy[:, j+1]
            end
            Phi_n = St + halpha2.*(a0[1:alpha_length, n+1].*cache.fy[:, 1] + zn_corr[:, n+1] + Phi)
            yn0 = y_pred
            fn0 = f_pred
            stop = false
            mu_it = 0

            yn1 = zeros(T, alpha_length)
            fn1 = zeros(T, length(yn1))
            while stop == false
                yn1 = Phi_n + halpha2.*fn0
                mu_it = mu_it + 1
                if mu == Inf
                    stop = norm(yn1-yn0, Inf) < abstol
                    if mu_it > 100 && ~stop
                        stop = 1
                    end
                else
                    stop = (mu_it == mu)
                end
                prob.f(fn1, yn1, p, mesh[n+1])
                yn0 = yn1; fn0 = fn1
            end
            cache.y[:, n+1] = yn1
            cache.fy[:, n+1] = fn1
        end
    end
end

function  starting_term(t, u0, m_alpha, t0, m_alpha_factorial)
    ys = zeros(size(u0, 1), 1)
    for k = 1 : maximum(m_alpha)
        if length(m_alpha) == 1
            ys = ys + (t-t0)^(k-1)/m_alpha_factorial[k]*u0[:, k]
        else
            i_alpha = findall(x -> x>=k, m_alpha)
            ys[i_alpha, 1] = ys[i_alpha, 1] + (t-t0)^(k-1)*u0[i_alpha, k]./m_alpha_factorial[i_alpha, k]
        end
    end
    return ys
end