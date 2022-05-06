"""
    solve(prob::MultiTermsFODEProblem, h, PIIMRect())

Use implicit product integration rectangular type method to solve multi-terms FODE.
"""
struct PIIMRect <: FractionalDiffEqAlgorithm end

function solve(prob::MultiTermsFODEProblem, h, ::PIIMRect)
    @unpack parameters, orders, rightfun, u0, t0, T = prob
    
    J_fun(x, y) = ForwardDiff.derivative(x -> rightfun(x, y), x)
    
    Q = length(orders)
    orders= sort(orders)
    i_al = sortperm(orders, rev=true)
    parameters = parameters[i_al]
    al_Q = orders[end]
    al_i = orders[1:end-1]
    lam_Q = parameters[end]
    lam_rat_i = parameters[1:end-1]/lam_Q
    m_Q = ceil(Int64, al_Q)
    m_i = ceil.(Int64, orders[1:end-1])
    bet = [al_Q .- al_i; al_Q]
    
    itmax = 100
    tol = 1.0e-6 
    
    gamma_val = zeros(Q, m_Q)
    for i = 1 : Q-1
        k = collect(Int, 0:m_i[i]-1)
        gamma_val[i, k.+1] = gamma.(k.+bet[i].+1)
    end
    k = collect(0:m_Q-1)
    gamma_val[Q, :] = factorial.(k)
    
    
    problem_size = size(u0, 1)
    
    r::Int = 64
    N::Int = ceil(Int, (T-t0)/h)
    Nr::Int = ceil(Int, (N+1)/r)*r
    Qr::Int = ceil(Int, log2((Nr)/r)) - 1
    NNr::Int = 2^(Qr+1)*r
    
    y = zeros(problem_size, N+1)
    fy = zeros(problem_size, N+1)
    zn = zeros(problem_size, NNr+1, Q)
    
    nvett = collect(0:NNr+1)
    bn = zeros(Q, NNr+1)
    for i = 1 : Q
        nbeta = nvett.^bet[i]
        bn[i, :] = (nbeta[2:end] - nbeta[1:end-1]) * h^bet[i]/gamma(bet[i]+1)
    end
    C = 0
    for i = 1 : Q-1
        C = C + lam_rat_i[i]*bn[i, 1]
    end
    
    t = collect(0:N)*h
    y[:, 1] = u0[:, 1]
    fy[:, 1] .= f_vectorfield(t0, u0[:, 1], rightfun)
    (y, fy) = PIExTriangolo(1, r-1, t, y, fy, zn, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, C, tol, itmax, J_fun)
    
    ff = zeros(1, 2^(Qr+2)); ff[1:2] = [0 2]; card_ff = 2
    nx0 = 0; nu0 = 0
    for qr = 0 : Qr
        L = 2^qr
        (y, fy) = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, nu0, t, y, fy, zn, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, C, J_fun, itmax, tol)
        ff[1:2*card_ff] = [ff[1:card_ff] ff[1:card_ff]]
        card_ff = 2*card_ff
        ff[card_ff] = 4*L
    end

    if T<t[N+1]
        c = (T - t[N])/h
        t[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end

    t = t[1:N+1] ; y = y[:, 1:N+1]
    return FODESolution(t, y[:])
end
    

function DisegnaBlocchi(L, ff, r, Nr, nx0, nu0, t, y, fy, zn, N , bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, C, J_fun, itmax, tol)
    
    nxi::Int = nx0
    nxf::Int = nx0 + L*r - 1
    nyi::Int = nu0
    nyf::Int = nu0 + L*r - 1
    is = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[1] = nxi
    s_nxf[1] = nxf
    s_nyi[1] = nyi
    s_nyf[1] = nyf
    
    i_triangolo = 0; stop = false
    while stop == false
        stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1 >= Nr-1)
        
        zn = PIExQuadrato(nxi, nxf, nyi, nyf, y, fy, zn, bn, problem_size, Q)
        
        (y, fy) = PIExTriangolo(nxi, nxi+r-1, t, y, fy, zn, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, C, tol, itmax, J_fun)
        i_triangolo = i_triangolo + 1
        
        if stop==false
            if nxi+r-1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1; nxf = s_nxf[is]  + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf
            else
                nxi = nxi + r; nxf = nxi + r - 1; nyi = nyf + 1; nyf = nyf + r
                is = is + 1
                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf
            end
        end
        
    end
    return y, fy
end

function PIExQuadrato(nxi, nxf, nyi, nyf, y, fy, zn, bn,  problem_size, Q)
    coef_beg = nxi-nyf; coef_end = nxf-nyi+1
    funz_beg = nyi+1; funz_end = nyf+1
    
    for i = 1:Q
        vett_coef = bn[i, coef_beg:coef_end]
        if nyi == 0
            if i < Q
                vett_funz = [zeros(problem_size, 1)  y[:, funz_beg+1:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
            else
                vett_funz = [zeros(problem_size, 1)  fy[:, funz_beg+1:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
            end
        else
            if i < Q
                vett_funz = [y[:, funz_beg:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
            else
                vett_funz = [fy[:, funz_beg:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
            end
        end
        zzn = real.(FastConv(vett_coef, vett_funz))
        zn[:, nxi+1:nxf+1, i] = zn[:, nxi+1:nxf+1, i] + zzn[:, nxf-nyf+1:end]
    end
    return zn
end
    
    
    
 
function PIExTriangolo(nxi, nxf, t, y, fy, zn, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, C, tol, itmax, J_fun)

    for n = nxi:min(N, nxf)
        St = PIExStartingTerm_Multi(t[n+1], t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
        Phi_n = copy(St)

        for i = 1:Q-1
            temp = zn[:, n+1, i]
            for j = nxi:n-1
                temp = temp + bn[i, n-j+1]*y[:, j+1]
            end
            Phi_n = Phi_n - lam_rat_i[i]*temp
        end
        temp = zn[:, n+1, Q]
        for j = nxi : n-1
            temp = temp + bn[Q, n-j+1]*fy[:, j+1]
        end
        Phi_n = Phi_n + temp/lam_Q


        yn0 = y[:, n]; fn0 = f_vectorfield(t[n+1], yn0, rightfun)
        Jfn0 = Jf_vectorfield(t[n+1], yn0, J_fun)
        Gn0 = (1+C)*yn0 .- bn[Q, 1]./lam_Q*fn0 .- Phi_n
        stop = false; it = 0
        
        while ~stop
            
            JGn0 = (1+C)*zeros(problem_size, problem_size)+I .- bn[Q, 1]/lam_Q*Jfn0
            global yn1 = yn0 - JGn0\Gn0
            global fn1 = f_vectorfield(t[n+1], yn1, rightfun)
            Gn1 = (1+C)*yn1 .- bn[Q, 1]/lam_Q*fn1 .- Phi_n
            it = it + 1
            
            stop = (norm(yn1-yn0, Inf) < tol) || (norm(Gn1, Inf) < tol)
            if it > itmax && ~stop
                @warn "Non Convergence"
                stop = true
            end
            
            yn0 = yn1 ; Gn0 = Gn1
            if ~stop
                Jfn0 = Jf_vectorfield(t[n+1], yn0, J_fun)
            end
            
        end
    
        y[:, n+1] = yn1
        fy[:, n+1] .= fn1
            
    end
    return y, fy
end


function FastConv(x, y)
    Lx = length(x); Ly = size(y, 2); problem_size = size(y, 1)

    r = Lx
    z = zeros(Number, problem_size, r)
    X = ourfft(x, r)
    for i = 1:problem_size
        Y = ourfft(y[i, :]', r)
        Z = X.*Y
        z[i, :] = ourifft(Z, r)
    end
    return z
end

function ourfft(x, n)
    s=length(x)
    x=x[:]
    if s > n
        return fft(x[1:n])
    elseif s < n
        return fft([x; zeros(n-s)])
    else
        return fft(x)
    end
end

function ourifft(x, n)
    s=length(x)
    x=x[:]
    if s > n
        return ifft(x[1:n])
    elseif s < n
        return ifft([x; zeros(n-s)])
    else
        return ifft(x)
    end
end

f_vectorfield(t, y, rightfun)=rightfun(t, y)
Jf_vectorfield(t, y, fun)=fun(t, y)

function PIExStartingTerm_Multi(t,t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
    ys = zeros(problem_size)

    for k = 0:m_Q-1
        ys = ys .+ (t-t0)^k./gamma_val[Q, k+1]*u0[:, k+1]
    end
    for i = 1 : Q-1
        temp = zeros(problem_size)
        for k = 0 : m_i[i]-1
            temp = temp .+ (t-t0)^(k+bet[i])/gamma_val[i, k+1]*u0[:, k+1]
        end
        ys = ys + lam_rat_i[i]*temp
    end
    return ys'
end