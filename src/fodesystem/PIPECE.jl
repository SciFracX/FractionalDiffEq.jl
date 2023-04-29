#=
"""
    solve(prob::FODESystem, h, PECE())

Use the Adams-Bashforth-Moulton method to solve the system of FODEs.

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
=#

mutable struct M
    an
    bn
    a0
    halpha1
    halpha2
    mu
    mu_tol
    r
    index_fft
    an_fft
    bn_fft
end
#TODO: Rename as PIPECE
function solve(prob::FODESystem, h, ::PECE)
    @unpack f, α, u0, tspan, p = prob
    t0 = tspan[1]; T = tspan[2]
    METH = M(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)# Initialization
    mu_tol = 1.0e-6
    mu = 1
    α = α[:]

    # issue [#64](https://github.com/SciFracX/FractionalDiffEq.jl/issues/64)
    max_order = findmax(α)[1]
    if max_order > 1
        @error "This method doesn't support high order FDEs"
    end


    # Check compatibility size of the problem with number of fractional orders
    alpha_length = length(α)
    problem_size = size(u0, 1)

    m_alpha = ceil.(Int, α)
    m_alpha_factorial = zeros(alpha_length, maximum(m_alpha))
    for i = 1 : alpha_length
        for j = 0 : m_alpha[i]-1
            m_alpha_factorial[i, j+1] = factorial(j)
        end
    end

    f_temp = zeros(size(u0[:, 1]))
    #f_temp = sysf_vectorfield(t0, u0[:, 1], f)
    f(f_temp, u0[:, 1], p, t0)

    r::Int = 16
    N::Int = ceil(Int64, (T-t0)/h)
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
        if α[i_alpha] == α[1:i_alpha-1]
            push!(find_alpha, i_alpha)
        end

        if isempty(find_alpha) == false
            bn[i_alpha, :] = bn[find_alpha[1], :]
            an[i_alpha, :] = an[find_alpha[1], :]
            a0[i_alpha, :] = a0[find_alpha[1], :]
        else
            nalpha = nvett.^α[i_alpha]
            nalpha1 = nalpha.*nvett
            bn[i_alpha, :] = nalpha[2:end] - nalpha[1:end-1]
            an[i_alpha, :] = [1; (nalpha1[1:end-2] - 2*nalpha1[2:end-1] + nalpha1[3:end]) ]
            a0[i_alpha, :] = [0; nalpha1[1:end-2]-nalpha[2:end-1].*(nvett[2:end-1].-α[i_alpha].-1)]
        end
    end
    METH.bn = bn; METH.an = an; METH.a0 = a0
    METH.halpha1 = h.^α./gamma.(α.+1)
    METH.halpha2 = h.^α./gamma.(α.+2)
    METH.mu = mu; METH.mu_tol = mu_tol

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
        
        bn_fft = zeros(Complex, alpha_length, index_fft[2, Qr+1]); an_fft = copy(bn_fft)
        for l = 1 : Qr+1
            coef_end = 2^l*r
            for i_alpha = 1 : alpha_length
                find_alpha = Float64[]
                if α[i_alpha] == α[1:i_alpha-1]
                    push!(find_alpha, i_alpha)
                end
                if isempty(find_alpha) == false
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = bn_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = an_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
                else
                    bn_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(METH.bn[i_alpha, 1:coef_end], coef_end)
                    an_fft[i_alpha, index_fft[1, l]:index_fft[2, l]] = ourfft(METH.an[i_alpha, 1:coef_end], coef_end)
                end
            end
        end
        METH.bn_fft = bn_fft ; METH.an_fft = an_fft ; METH.index_fft = index_fft ;
    end

    # Initializing solution and proces of computation
    t = t0 .+ collect(0:N)*h
    y[:, 1] = u0[:, 1]
    fy[:, 1] = f_temp
    (y, fy) = ABM_triangolo(1, r-1, t, y, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, p, f) ;

    # Main process of computation by means of the FFT algorithm
    ff = zeros(1, 2^(Qr+2)); ff[1:2] = [0; 2] ; card_ff = 2
    nx0::Int = 0; ny0::Int = 0
    for qr = 0 : Qr
        L = 2^qr 
        (y, fy) = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, p, f)
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


function DisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, p, f)

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
        
        (zn_pred, zn_corr) = ABM_quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length)
        
        (y, fy) = ABM_triangolo(nxi, nxi+r-1, t, y, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, p, f)
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
    return y, fy
end

function ABM_quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length)
    coef_end::Int = nxf-nyi+1
    i_fft::Int = log2(coef_end/METH.r) 
    funz_beg::Int = nyi+1
    funz_end::Int = nyf+1
    Nnxf::Int = min(N, nxf)

    # Evaluation convolution segment for the predictor
    vett_funz = fy[:, funz_beg:funz_end]
    vett_funz_fft = rowfft(vett_funz, coef_end)
    zzn_pred = zeros(problem_size, coef_end)
    for i = 1 : problem_size
        i_alpha = min(alpha_length,i)
        Z = METH.bn_fft[i_alpha, METH.index_fft[1, i_fft]:METH.index_fft[2, i_fft]].*vett_funz_fft[i, :]
        zzn_pred[i, :] = real.(ourifft(Z, coef_end))
    end
    zzn_pred = zzn_pred[:, nxf-nyf:end-1]
    zn_pred[:, nxi+1:Nnxf+1] = zn_pred[:, nxi+1:Nnxf+1] + zzn_pred[:, 1:Nnxf-nxi+1]

    # Evaluation convolution segment for the corrector
    if METH.mu > 0
        if nyi == 0 # Evaluation of the lowest square
            vett_funz = [zeros(problem_size, 1) fy[:, funz_beg+1:funz_end]]
            vett_funz_fft = rowfft(vett_funz, coef_end)
        end
        zzn_corr = zeros(problem_size, coef_end)
        for i = 1 : problem_size
            i_alpha = min(alpha_length,i)
            Z = METH.an_fft[i_alpha,METH.index_fft[1, i_fft]:METH.index_fft[2, i_fft]].*vett_funz_fft[i, :]
            zzn_corr[i, :] = real.(ourifft(Z, coef_end))
        end
        zzn_corr = zzn_corr[:, nxf-nyf+1:end]
        zn_corr[:, nxi+1:Nnxf+1] = zn_corr[:, nxi+1:Nnxf+1] + zzn_corr[:, 1:Nnxf-nxi+1]
    else
        zn_corr = 0
    end
    return zn_pred, zn_corr
end



function ABM_triangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, METH, problem_size, alpha_length, m_alpha, m_alpha_factorial, u0, t0, p, f)

    for n = nxi:min(N, nxf)
        
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
        St = starting_term(t[n+1], u0, m_alpha, t0, m_alpha_factorial)
        y_pred = St + METH.halpha1.*(zn_pred[:, n+1] + Phi)
        f_pred = zeros(length(y_pred))
        #f_pred = sysf_vectorfield(t[n+1], y_pred, f)
        f(f_pred, y_pred, p, t[n+1])
        
        # Evaluation of the corrector
        if METH.mu == 0
            y[:, n+1] = y_pred
            fy[:, n+1] = f_pred
        else
            j_beg = nxi
            Phi = zeros(problem_size, 1)
            for j = j_beg : n-1
                Phi = Phi + METH.an[1:alpha_length, n-j+1].*fy[:, j+1]
            end
            Phi_n = St + METH.halpha2.*(METH.a0[1:alpha_length, n+1].*fy[:, 1] + zn_corr[:, n+1] + Phi)
            yn0 = y_pred
            fn0 = f_pred
            stop = false
            mu_it = 0
            while stop == false
                global yn1 = Phi_n + METH.halpha2.*fn0
                mu_it = mu_it + 1
                if METH.mu == Inf
                    stop = norm(yn1-yn0, Inf) < METH.mu_tol
                    if mu_it > 100 && ~stop
                        stop = 1
                    end
                else
                    stop = (mu_it == METH.mu)
                end
                global fn1 = zeros(length(yn1))
                #sysf_vectorfield(t[n+1], yn1, f)
                f(fn1, yn1, p, t[n+1])
                yn0 = yn1; fn0 = fn1
            end
            y[:, n+1] = yn1
            fy[:, n+1] = fn1
        end
    end
    return y, fy
end


sysf_vectorfield(t, y, f_fun) = f_fun(t, y)

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