@concrete mutable struct PITrapCache{iip, T}
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
    Jfdefun

    r
    N
    Nr
    Qr
    NNr

    an
    a0

    halpha2
    abstol
    maxiters
    index_fft
    an_fft

    kwargs
end

Base.eltype(::PITrapCache{iip, T}) where {iip, T} = T

struct PITrap <: FODEAlgorithm end

function SciMLBase.__init(prob::FODEProblem, alg::PITrap; dt = 0.0, abstol = 1e-6, maxiters = 1000, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    prob = _is_need_convert!(prob)
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    order = order[:]

    alpha_length = length(order)
    problem_size = length(u0)

    m_alpha = ceil.(Int, order)
    m_alpha_factorial = zeros(alpha_length, maximum(m_alpha))
    for i = 1 : alpha_length
        for j = 0 : m_alpha[i]-1
            m_alpha_factorial[i, j+1] = factorial(j)
        end
    end

    f_temp = zeros(length(u0[:, 1]))
    f(f_temp, u0[:, 1], p, t0)

    r = 16
    N = ceil(Int64, (tfinal-t0)/dt)
    Nr = ceil(Int64, (N+1)/r)*r
    Qr = ceil(Int64, log2(Nr/r)) - 1
    NNr = 2^(Qr+1)*r

    # Preallocation of some variables
    y = zeros(T, problem_size, N+1)
    fy = zeros(T, problem_size, N+1)
    zn = zeros(problem_size, NNr+1)

    # generate jacobian of the input function
    Jfdefun(t, u) = jacobian_of_fdefun(prob.f, t, u, p)

    # Evaluation of coefficients of the PECE method
    nvett = 0:NNr+1 |> collect
    an = zeros(alpha_length, NNr+1)
    a0 = copy(an)
    for i_alpha = 1:alpha_length
        find_alpha = Float64[]
        if order[i_alpha] == order[1:i_alpha-1]
            push!(find_alpha, i_alpha)
        end

        if isempty(find_alpha) == false
            an[i_alpha, :] = an[find_alpha[1], :]
            a0[i_alpha, :] = a0[find_alpha[1], :]
        else
            nalpha = nvett.^order[i_alpha]
            nalpha1 = nalpha.*nvett
            an[i_alpha, :] = [ 1; (nalpha1[1:end-2] - 2*nalpha1[2:end-1] + nalpha1[3:end]) ] ;
            a0[i_alpha, :] = [ 0; nalpha1[1:end-2]-nalpha[2:end-1].*(nvett[2:end-1] .- order[i_alpha] .-1)]
        end
    end
    halpha2 = dt.^order./gamma.(order.+2)

    
    if Qr >= 0
        index_fft = zeros(Int, 2, Qr+1)
        for l = 1 : Qr+1
            if l == 1
                index_fft[1, l] = 1; index_fft[2, l] = r*2
            else
                index_fft[1, l] = index_fft[2, l-1]+1; index_fft[2, l] = index_fft[2, l-1]+2^l*r
            end
        end
        
        an_fft = zeros(Complex, alpha_length, index_fft[2, Qr+1])
        for l = 1:Qr+1
            coef_end = 2^l*r
            for i_alpha = 1 : alpha_length
                find_alpha = Float64[]
                if order[i_alpha] == order[1:i_alpha-1]
                    push!(find_alpha, i_alpha)
                end
                if isempty(find_alpha) == false
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = an_fft[find_alpha[1], index_fft[1, l]:index_fft[2, l]]
                else
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(an[i_alpha, 1:coef_end], coef_end)
                end
            end
        end
    else
        index_fft = 0
        an_fft = 0
    end
    # Initializing solution and proces of computation
    mesh = t0 .+ collect(0:N)*dt
    y[:, 1] = u0[:, 1]
    fy[:, 1] = f_temp
    return PITrapCache{iip, T}(prob, alg, mesh, u0, order, m_alpha, m_alpha_factorial, y, fy, p, problem_size,
                             zn, Jfdefun, r, N, Nr, Qr, NNr, an, a0, halpha2, abstol, maxiters, index_fft, an_fft, kwargs)
end

function SciMLBase.solve!(cache::PITrapCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, y, fy, r, N, Nr, Qr, NNr, an, a0, halpha2, abstol, index_fft, an_fft, kwargs = cache
    t0 = mesh[1]
    tfinal = mesh[end]
    PITrap_triangolo(cache, 1, r-1)

    # Main process of computation by means of the FFT algorithm
    ff = zeros(1, 2^(Qr+2)); ff[1:2] = [0; 2] ; card_ff = 2
    nx0 = 0; ny0 = 0
    for qr = 0 : Qr
        L = 2^qr 
        PITrap_disegna_blocchi(cache, L, ff, nx0+L*r, ny0)
        ff[1:2*card_ff] = [ff[1:card_ff]; ff[1:card_ff]] 
        card_ff = 2*card_ff
        ff[card_ff] = 4*L
    end

    # Evaluation solution in tfinal when tfinal is not in the mesh
    if tfinal < mesh[N+1]
        c = (tfinal - mesh[N])/dt
        mesh[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end

    mesh = mesh[1:N+1]; y = cache.y[:, 1:N+1]
    u = collect(Vector{eltype(u0)}, eachcol(y))

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end


function PITrap_disegna_blocchi(cache::PITrapCache{iip, T}, L::P, ff, nx0::P, ny0::P) where {P <: Integer, iip, T}
    @unpack mesh, N, r, Nr = cache

    nxi::Int = nx0; nxf::Int = nx0 + L*r - 1
    nyi::Int = ny0; nyf::Int = ny0 + L*r - 1
    is = 1
    s_nxi = zeros(T, N)
    s_nxf = zeros(T, N)
    s_nyi = zeros(T, N)
    s_nyf = zeros(T, N)
    s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf

    i_triangolo = 0; stop = false
    while stop == false
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1)
        
        PITrap_quadrato(cache, nxi, nxf, nyi, nyf)
        PITrap_triangolo(cache, nxi, nxi+r-1)
        i_triangolo = i_triangolo + 1
        
        if stop == false
            if nxi+r-1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is]  + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            else
                nxi = nxi + r; nxf = nxi + r - 1; nyi = nyf + 1; nyf = nyf + r
                is = is + 1
                s_nxi[is] = nxi; s_nxf[is] = nxf; s_nyi[is] = nyi; s_nyf[is] = nyf
            end
        end
    end
end

function PITrap_quadrato(cache::PITrapCache{iip, T}, nxi::P, nxf::P, nyi::P, nyf::P) where {P <: Integer, iip, T}
    @unpack prob, mesh, r, N, Nr, Qr, NNr, problem_size, an, a0, halpha2, abstol, index_fft, an_fft = cache
    coef_end = nxf-nyi+1
    alpha_length = length(prob.order)
    i_fft::Int = log2(coef_end/r) 
    funz_beg = nyi+1; funz_end = nyf+1
    Nnxf = min(N, nxf)

    if nyi == 0
        vett_funz = [zeros(problem_size, 1) cache.fy[:, funz_beg+1:funz_end]]
    else
        vett_funz = cache.fy[:, funz_beg:funz_end]
    end
    vett_funz_fft = rowfft(vett_funz, coef_end)
    zzn = zeros(problem_size, coef_end)
    for i = 1 : problem_size
        i_alpha::Int = min(alpha_length, i)
        if abs(prob.order[i_alpha]-1)>1e-14
        Z = an_fft[i_alpha, index_fft[1, i_fft]:index_fft[2, i_fft]].*vett_funz_fft[i, :]
        zzn[i, :] = real.(ourifft(Z, coef_end))
        end
    end
    zzn = zzn[:, nxf-nyf+1:end]
    cache.zn[:, nxi+1:Nnxf+1] = cache.zn[:, nxi+1:Nnxf+1] + zzn[:, 1:Nnxf-nxi+1]
end



function PITrap_triangolo(cache::PITrapCache{iip, T}, nxi::P, nxf::P) where {P <: Integer, iip, T}
    @unpack prob, mesh, u0, order, m_alpha, m_alpha_factorial, p, problem_size, zn, Jfdefun, N, an, a0, halpha2, abstol, maxiters, index_fft, an_fft = cache

    alpha_length = length(order)
    for n = nxi:min(N, nxf)
        n1 = n+1
        St = PITrap_system_starting_term(cache, mesh[n+1])
        # Evaluation of the predictor
        Phi = zeros(problem_size, 1)
        for j = nxi:n-1
            Phi = Phi + an[1:alpha_length,n-j+1].*cache.fy[:, j+1]
        end
        Phi_n = St + halpha2 .*(a0[1:alpha_length, n+1] .* cache.fy[:, 1] + cache.zn[:, n+1] + Phi)
        
        yn0 = cache.y[:, n]
        fn0 = zeros(T, problem_size); Jfn0 = zeros(T, problem_size, problem_size)
        prob.f(fn0, yn0, p, mesh[n+1])
        Jfn0 = Jf_vectorfield(mesh[n+1], yn0, Jfdefun)
        Gn0 = yn0 - halpha2 .*fn0 - Phi_n

        stop = false; it = 0
        yn1 = similar(yn0)
        fn1 = similar(yn0)

        while ~stop
            JGn0 = zeros(T, problem_size, problem_size)+I - diagm(halpha2)*Jfn0
            yn1 = yn0 - JGn0\Gn0
            prob.f(fn1, yn1, p, mesh[n+1])
            Gn1 = yn1 - halpha2.*fn1 - Phi_n
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

function  PITrap_system_starting_term(cache::PITrapCache{iip, T}, t) where {iip, T}
    @unpack mesh, u0, m_alpha, m_alpha_factorial = cache
    t0 = mesh[1]
    ys = zeros(length(u0))
    for k = 1 : maximum(m_alpha)
        if length(m_alpha) == 1
            ys = ys .+ (t-t0)^(k-1)/m_alpha_factorial[k]*u0[k]
        else
            i_alpha = findall(x -> x>=k, m_alpha)
            ys[i_alpha] = ys[i_alpha] + (t-t0)^(k-1)*u0[i_alpha, k]./m_alpha_factorial[i_alpha, k]
        end
    end
    return ys
end