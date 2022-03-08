
function PIEX(al,lam,f_fun,t0,T,y0,h)
    Q = length(al)
    al= sort(al)
    i_al = sortperm(al)
    lam = lam[i_al]
    al_Q = al[end]
    al_i = al[1:end-1]
    lam_Q = lam[end]
    lam_rat_i = lam[1:end-1]/lam_Q
    m_Q::Int64 = ceil(al_Q)
    m_i = ceil.(al[1:end-1])
    bet = [al_Q .- al_i ; al_Q]
    
    gamma_val = zeros(Q, m_Q)
    for i = 1 : Q-1
        k = collect(Int, 0:m_i[i]-1)
        gamma_val[i, k.+1] = gamma.(k.+bet[i].+1)
    end
    k = collect(0:m_Q-1)
    gamma_val[Q, :] = factorial.(k)
    
    
    problem_size = size(y0, 1)
    
    r = 16
    N = ceil(Int, (T-t0)/h)
    Nr = ceil(Int, (N+1)/r)*r
    Qr = ceil(Int, log2((Nr)/r)) - 1
    NNr = 2^(Qr+1)*r
    
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
    y[:, 1] = y0[:, 1]
    fy[:, 1] .= f_vectorfield(t0, y0[:, 1], f_fun)
    (y, fy) = Triangolo(1, r-1, t, y, fy, zn, N, bn, t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, f_fun, lam_Q)
    
    ff = zeros(1, 2^(Qr+2)) ; ff[1:2] = [0 2] ; card_ff = 2
    nx0 = 0 ; ny0 = 0 ;
    for qr = 0 : Qr
        L = 2^qr ; 
        (y, fy) = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, zn, N, bn, t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, f_fun, lam_Q) ;
        ff[1:2*card_ff] = [ff[1:card_ff] ff[1:card_ff]]
        card_ff = 2*card_ff
        ff[card_ff] = 4*L
    end
    
    if T < t[N+1]
        c = (T - t[N])/h
        t[N+1] = tfinal
        y[:, N+1] = (1-c)*y[:, N] + c*y[:, N+1]
    end
    
    t = t[1:N+1] ; y = y[:, 1:N+1]
    return t, y
end
    

function DisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, zn, N , bn, t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, f_fun, lam_Q)
    
    nxi = nx0 ; nxf = nx0 + L*r - 1
    nyi = ny0 ; nyf = ny0 + L*r - 1
    is = 1
    s_nxi = zeros(N)
    s_nxf = zeros(N)
    s_nyi = zeros(N)
    s_nyf = zeros(N)
    s_nxi[is] = nxi
    s_nxf[is] = nxf
    s_nyi[is] = nyi
    s_nyf[is] = nyf
    
    i_triangolo = 0; stop = 0
    while stop!==0
        
        stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1)
        
        zn = Quadrato(nxi, nxf, nyi, nyf, y, fy, zn, bn, problem_size, Q) ;
        
        (y, fy) = Triangolo(nxi, nxi+r-1, t, y, fy, zn, N, bn, t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, f_fun, lam_Q) ;
        i_triangolo = i_triangolo + 1
        
        if ~stop
            if nxi+r-1 == nxf
                i_Delta = ff[i_triangolo]
                Delta = i_Delta*r
                nxi = s_nxf[is]+1 ; nxf = s_nxf[is]  + Delta
                nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf
            else
                nxi = nxi + r ; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf = nyf + r  ;
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
    
 
function Quadrato(nxi, nxf, nyi, nyf, y, fy, zn, bn,  problem_size, Q)
    
    coef_beg = nxi-nyf ; coef_end = nxf-nyi+1
    funz_beg = nyi+1 ; funz_end = nyf+1
    
    
    for i = 1:Q
        vett_coef = bn[i, coef_beg:coef_end]
        if i < Q
            vett_funz = [y[:, funz_beg:funz_end], zeros(problem_size, funz_end-funz_beg+1)]
        else
            vett_funz = [fy[:, funz_beg:funz_end], zeros(problem_size, funz_end-funz_beg+1)]
        end
        zzn = real(FastConv(vett_coef,vett_funz))
        zn[:, nxi+1:nxf+1, i] = zn[:, nxi+1:nxf+1, i] + zzn[:, nxf-nyf+1-1:end-1]
    end
    return zn
end
    
    
    
 
function Triangolo(nxi, nxf, t, y, fy, zn, N, bn, t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val, f_fun, lam_Q)

    for n = nxi : min(N, nxf)
        St = StartingTerm_Multi(t[n+1], t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
        Phi_n = St
        if nxi == 1
            j_beg = 0
        else
            j_beg = nxi
        end
        
        for i = 1 : Q-1
            temp = zn[:, n+1, i]
            for j = j_beg:n-1
                temp = temp + bn[i, n-j]*y[:, j+1]
            end
            Phi_n = Phi_n - lam_rat_i[i]*temp
        end
        temp = zn[:, n+1, Q]
        for j = j_beg : n-1
            temp = temp + bn[Q, n-j]*fy[:, j+1]
        end
        Phi_n = Phi_n + temp/lam_Q
        
        y[:, n+1] = Phi_n
        fy[:, n+1] .= f_vectorfield(t[n+1], y[:, n+1], f_fun) 
            
    end
    return y, fy
end
    
    
function FastConv(x, y)
    Lx = length(x); Ly = size(y, 2) ; problem_size = size(y, 1)

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
    
    
f_vectorfield(t,y,f_fun)=f_fun(t, y)
    
function StartingTerm_Multi(t,t0, problem_size, y0, Q, m_Q, m_i, bet, lam_rat_i, gamma_val)
    ys = zeros(problem_size, 1)

    for k = 0:m_Q-1
        ys = ys .+ (t-t0)^k./gamma_val[Q, Int64(k+1)]*y0[:, k+1]
    end
    for i = 1 : Q-1
        temp = zeros(problem_size, 1)
        for k = 0 : Int64(m_i[i]-1)
            temp = temp + (t-t0)^(k+bet[i])/gamma_val[i, k+1]*y0[:, k+1]
        end
        ys = ys + lam_rat_i[i]*temp
    end
    return ys
end
#=
test(t, y)=172/125*cos(4/5*t)
(t, y)=PIEX([3, 2.5, 2, 1, 0.5, 1], [1, 1/16, 4/5, 3/2, 1/25, 6/5], test, 0, 20, [0 0 0 0 0 0], 0.01)
=#