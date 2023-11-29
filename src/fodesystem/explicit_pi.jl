

"""
    solve(prob::FODEProblem, h, PIEX())

Use explicit Product integration method to solve system of FODE.
"""
function solve(prob::FODEProblem, h, ::PIEX)
    @unpack f, order, u0, tspan, p = prob

    t0 = tspan[1]; T = tspan[2]

    # issue [#64](https://github.com/SciFracX/FractionalDiffEq.jl/issues/64)
    max_order = findmax(order)[1]
    if max_order > 1
        @error "This method doesn't support high order FDEs"
    end

    METH = M(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)# Initialization
    order = order[:]
    # Check compatibility size of the problem with number of fractional orders
    alpha_length = length(order)
    problem_size = size(u0, 1)

    m_alpha = ceil.(Int, order)
    m_alpha_factorial = zeros(alpha_length, maximum(m_alpha))
    for i = 1 : alpha_length
        for j = 0 : m_alpha[i]-1
            m_alpha_factorial[i, j+1] = factorial(j)
        end
    end

    f_temp = zeros(length(u0[:, 1]))
    f(f_temp, u0[:, 1], p, t0)

    r::Int = 16
    N::Int = ceil(Int64, (T-t0)/h)
    Nr::Int = ceil(Int64, (N+1)/r)*r
    Qr::Int = ceil(Int64, log2(Nr/r)) - 1
    NNr::Int = 2^(Qr+1)*r

    # Preallocation of some variables
    y = zeros(problem_size, N+1)
    fy = zeros(problem_size, N+1)
    zn = zeros(problem_size, NNr+1)

    # Evaluation of coefficients of the PECE method
    nvett = 0:NNr+1 |> collect
    bn = zeros(alpha_length, NNr+1)#; an = copy(bn); a0 = copy(bn)
    for i_alpha = 1:alpha_length
        find_alpha = Float64[]
        if order[i_alpha] == order[1:i_alpha-1]
            push!(find_alpha, i_alpha)
        end

        if isempty(find_alpha) == false
            bn[i_alpha, :] = bn[find_alpha[1], :]
        elseif abs(order[i_alpha]-1) < 1e-14
            bn[i_alpha, :] = [1; zeros(NNr)]
        else
            nalpha = nvett.^order[i_alpha]
            bn[i_alpha, :] = nalpha[2:end] - nalpha[1:end-1]
        end
    end
    METH.bn = bn
    METH.halpha1 = h.^order./gamma.(order.+1)

    # Evaluation of FFT of coefficients of the PECE method
    METH.r = r 
    if Qr >= 0
        index_fft = zeros(Int, 2, Qr+1)
        for l = 1 : Qr+1
            if l == 1
                index_fft[1, l] = 1; index_fft[2, l] = r*2
            else
                index_fft[1, l] = index_fft[2, l-1]+1; index_fft[2, l] = index_fft[2, l-1]+2^l*r
            end
        end
        
        bn_fft = zeros(Complex, alpha_length, index_fft[2, end]);# an_fft = copy(bn_fft)
        for l = 1:Qr+1
            coef_end = 2^l*r
            for i_alpha = 1 : alpha_length
                find_alpha = Float64[]
                if order[i_alpha] == order[1:i_alpha-1]
                    push!(find_alpha, i_alpha)
                end
                if isempty(find_alpha) == false
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = bn_fft[find_alpha[1], index_fft[1, l]:index_fft[2, l]]
                else
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(METH.bn[i_alpha, 1:coef_end], coef_end)
                end
            end
        end
        METH.bn_fft = bn_fft ; METH.index_fft = index_fft
    end

    # Initializing solution and proces of computation
    t = t0 .+ collect(0:N)*h
    y[:, 1] = u0[:, 1]
    fy[:, 1] = f_temp
    (y, fy) = PIEX_system_triangolo(1, r-1, t, y, fy, zn, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, f, order, p)

    # Main process of computation by means of the FFT algorithm
    ff = zeros(1, 2^(Qr+2)); ff[1:2] = [0; 2] ; card_ff = 2
    nx0::Int = 0; ny0::Int = 0
    for qr = 0 : Qr
        L = 2^qr 
        (y, fy) = PIEX_system_disegna_blocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, zn, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, f, order, p)
        ff[1:2*card_ff] = [ff[1:card_ff]; ff[1:card_ff]] 
        card_ff = 2*card_ff
        ff[card_ff] = 4*L
    end

    # Evaluation solution in T when T is not in the mesh
    if T < t[N+1]
        c = (T - t[N])/h
        t[N+1] = T
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end

    t = t[1:N+1]; y = y[:, 1:N+1]
    return FODESystemSolution(t, y)
end


function PIEX_system_disegna_blocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, zn, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, f, alpha, p)

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
        
        zn = PIEX_system_quadrato(nxi, nxf, nyi, nyf, fy, zn, N, METH, problem_size, alpha_length, alpha)
        
        (y, fy) = PIEX_system_triangolo(nxi, nxi+r-1, t, y, fy, zn, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, f, alpha, p)
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
    return y, fy
end

function PIEX_system_quadrato(nxi, nxf, nyi, nyf, fy, zn, N, METH, problem_size, alpha_length, alpha)
    coef_end::Int = nxf-nyi+1
    i_fft::Int = log2(coef_end/METH.r) 
    funz_beg::Int = nyi+1; funz_end::Int = nyf+1
    Nnxf::Int = min(N, nxf)

    # Evaluation convolution segment for the predictor
    vett_funz = fy[:, funz_beg:funz_end]
    vett_funz_fft = rowfft(vett_funz, coef_end)
    zzn = zeros(problem_size, coef_end)
    for i = 1 : problem_size
        i_alpha = min(alpha_length, i)
        if abs(alpha[i_alpha]-1)>1e-14
        Z = METH.bn_fft[i_alpha, METH.index_fft[1, i_fft]:METH.index_fft[2, i_fft]].*vett_funz_fft[i, :]
        zzn[i, :] = real.(ourifft(Z, coef_end))
        end
    end
    zzn = zzn[:, nxf-nyf:end-1]
    zn[:, nxi+1:Nnxf+1] = zn[:, nxi+1:Nnxf+1] + zzn[:, 1:Nnxf-nxi+1]
    return zn
end



function PIEX_system_triangolo(nxi, nxf, t, y, fy, zn, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, f, alpha, p)

    for n = nxi:min(N, nxf)
        St = PIEX_system_starting_term(t[n+1], u0, m_alpha, t0, m_alpha_factorial)        
        # Evaluation of the predictor
        Phi = zeros(problem_size, 1)
        if nxi == 1 # Case of the first triangle
            j_beg::Int = 0
        else # Case of any triangle but not the first
            j_beg = nxi
        end
        for j = j_beg:n-1
            Phi = Phi + METH.bn[1:alpha_length,n-j].*fy[:, j+1]
        end

        i_alpha_1 = findall(alpha -> abs(alpha - 1) < 1e-14, alpha)
        Phi[i_alpha_1] = fy[i_alpha_1, n]
        St[i_alpha_1] = y[i_alpha_1, n]
        
        y[:, n+1] = St + METH.halpha1.*(zn[:, n+1] + Phi)
        temp = zeros(length(y[:, n+1]))
        f(temp, y[:, n+1], p, t[n+1])
        fy[:, n+1] = temp
        
    end
    return y, fy
end

function  PIEX_system_starting_term(t, u0, m_alpha, t0, m_alpha_factorial)
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