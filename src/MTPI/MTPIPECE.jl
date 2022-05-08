"""
    solve(prob::MultiTermsFODEProblem, h, PIPECE())

Use product integration predictor-corrector method to solve multi-terms FODE.
"""
struct PIPECE <: FractionalDiffEqAlgorithm end

function solve(prob::MultiTermsFODEProblem, h, ::PIPECE)
    @unpack parameters, orders, rightfun, u0, t0, T = prob
    mu=1
    mu_tol=1e-6
    u0 = u0[:]'
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
    zn_pred = zeros(problem_size, NNr+1, Q)
    if mu>0
        zn_corr = zeros(problem_size, NNr+1, Q)
    else
        zn_corr = 0
    end
    
    nvett = collect(0:NNr+1)
    bn = zeros(Q, NNr+1); an = zeros(Q, NNr+1); a0 = zeros(Q, NNr+1)
    for i = 1 : Q
        nbeta = nvett.^bet[i]; nbeta1 = nbeta.*nvett
        bn[i, :] = (nbeta[2:end] - nbeta[1:end-1]) * h^bet[i] / gamma(bet[i]+1)
        an[i, :] = [1;  (nbeta1[1:end-2] - 2*nbeta1[2:end-1] + nbeta1[3:end]) ]*h^bet[i] / gamma(bet[i]+2)
        a0[i, :] = [0;  nbeta1[1:end-2] - nbeta[2:end-1].*(nvett[2:end-1].-bet[i].-1)]*h^bet[i] / gamma(bet[i]+2)
    end
    C = 0
    for i = 1 : Q-1
        C = C + lam_rat_i[i]*an[i, 1]
    end
    
    t = collect(0:N)*h
    y[:, 1] = u0[:, 1]
    fy[:, 1] .= f_vectorfield(t0, u0[:, 1], rightfun)
    (y, fy) = PIPECETriangolo(1, r-1, t, y, fy, zn_pred, zn_corr, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, an, a0, mu, mu_tol, C)
    
    ff = zeros(1, 2^(Qr+2)); ff[1:2] = [0 2]; card_ff = 2
    nx0 = 0; nu0 = 0
    for qr = 0 : Qr
        L = 2^qr
        (y, fy) = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, nu0, t, y, fy, zn_pred, zn_corr, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, an, a0, mu, mu_tol, C)
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
    

function DisegnaBlocchi(L, ff, r, Nr, nx0, nu0, t, y, fy, zn_pred, zn_corr, N , bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, an, a0, mu, mu_tol, C)
    
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
        
        (zn_pred, zn_corr) = PIPECEQuadrato(nxi, nxf, nyi, nyf, y, fy, zn_pred, zn_corr, bn, problem_size, Q, an, mu)
        
        (y, fy) = PIPECETriangolo(nxi, nxi+r-1, t, y, fy, zn_pred, zn_corr, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, an, a0, mu, mu_tol, C)
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

function PIPECEQuadrato(nxi, nxf, nyi, nyf, y, fy, zn_pred, zn_corr, bn,  problem_size, Q, an, mu)
    coef_beg = nxi-nyf; coef_end = nxf-nyi+1
    funz_beg = nyi+1; funz_end = nyf+1
    
    for i = 1:Q
        vett_coef = bn[i, coef_beg:coef_end]
        if i < Q
            vett_funz = [y[:, funz_beg:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
        else
            vett_funz = [fy[:, funz_beg:funz_end] zeros(problem_size, funz_end-funz_beg+1)]
        end
        zzn_pred = real.(FastConv(vett_coef, vett_funz))
        zn_pred[:, nxi+1:nxf+1, i] = zn_pred[:, nxi+1:nxf+1, i] + zzn_pred[:, nxf-nyf:end-1]
    end

    if mu > 0
        for i = 1 : Q
            vett_coef = an[i,coef_beg:coef_end]
            if nyi == 0
                if i < Q
                    vett_funz = [zeros(problem_size, 1)  y[:, funz_beg+1:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
                else
                    vett_funz = [zeros(problem_size, 1)  fy[:, funz_beg+1:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
                end
            else
                if i < Q
                    vett_funz = [y[:, funz_beg:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
                else
                    vett_funz = [fy[:, funz_beg:funz_end]  zeros(problem_size, funz_end-funz_beg+1)]
                end
            end
            zzn_corr = real(FastConv(vett_coef, vett_funz))
            zn_corr[:, nxi+1:nxf+1, i] = zn_corr[:, nxi+1:nxf+1, i] + zzn_corr[:, nxf-nyf+1:end]
        end
    else
        zn_corr = 0
    end

    return zn_pred, zn_corr
end
    
    
    
 
function PIPECETriangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, bn, t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, rightfun, lam_Q, an, a0, mu, mu_tol, C)

    for n = nxi:min(N, nxf)
        St = PIPECEStartingTerm_Multi(t[n+1], t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
        
        Phi_n = copy(St)
        if nxi == 1
            j_beg = 0
        else
            j_beg = nxi
        end
        
        for i = 1:Q-1
            temp = zn_pred[:, n+1, i]
            for j = j_beg:n-1
                temp = temp + bn[i, n-j]*y[:, j+1]
            end
            Phi_n = Phi_n - lam_rat_i[i]*temp
        end
        temp = zn_pred[:, n+1, Q]
        for j = j_beg : n-1
            temp = temp + bn[Q, n-j]*fy[:, j+1]
        end
        Phi_n = Phi_n + temp/lam_Q
        y_pred = copy(Phi_n)
        f_pred = copy(f_vectorfield(t[n+1], y_pred, rightfun))

        if mu == 0
            y[:,n+1] = y_pred
            fy[:, n+1] = f_pred
        else
            j_beg = nxi
            Phi_n = copy(St)
            for i = 1 : Q-1
                temp = a0[i, n+1]*y[:, 1] + zn_corr[:,n+1,i]
                for j = j_beg : n-1
                    temp += an[i,n-j+1]*y[:,j+1]#Possiable bugs
                end
                Phi_n = Phi_n - lam_rat_i[i]*temp
            end
            temp = a0[Q, n+1]*fy[:, 1] + zn_corr[:, n+1, Q]
            for j = j_beg : n-1
                temp = temp + an[Q, n-j+1]*fy[:, j+1]
            end
            Phi_n = Phi_n + temp/lam_Q
            
            yn0 = y_pred ; fn0 = f_pred
            stop = false; mu_it = 0
            while ~stop
                global yn1 = Phi_n .- C*yn0 .+ an[Q, 1]./lam_Q*fn0
                mu_it = mu_it + 1
                if mu == Inf
                    stop = (norm(yn1-yn0,inf) < mu_tol)
                    if mu_it > 100 && ~stop
                        @warn "Non convergence"
                        stop = 1
                    end
                else
                    stop = mu_it == mu
                end
                global fn1 = f_vectorfield(t[n+1], yn1, rightfun)           
                yn0 = yn1; fn0 = fn1
            end

        y[:, n+1] = yn1
        fy[:, n+1] .= fn1
        end
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
function process_rightfun(t, y, rightfun)
    if typeof(rightfun) <: Function
        return rightfun(t, y)
    else
        return rightfun*ones(length(y))
    end
end

function PIPECEStartingTerm_Multi(t,t0, problem_size, u0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
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