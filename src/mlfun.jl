## Classical and two-parametric Mittag Leffler function algorithms is in
## Rudolfo Gorenflo, Joulia Loutchko and Yuri Loutchko,
## *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal, **(2002)**
"""
```tex
@ARTICLE{Gorenflo_computationof,
    author = {Rudolf Gorenflo and Joulia Loutchko and Yuri Luchko},
    title = {Computation of the Mittag-Leffler function Eα, β(z) and its derivative},
    journal = {FRACT. CALC. APPL. ANAL},
    year = {},
    pages = {2002}
}
```
"""
# Classical and Two-Parametric Mittag Leffler function is modified from [John Lapeyre](https://github.com/jlapeyre)'s [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl) module.
# Since Mittag Leffler functions is widely used in Fractional Differential Equations, so we decided to has the Mittag Leffler function build in.

function ourquadgk(f, a, b)
    (result, _) = quadgk(f, a, b; order=17)
    return result
end

function P(α, β, ϵ, ϕ, z)
    ω = ϕ*(1 + (1 - β)/α) + ϵ^(1/α)*sin(ϕ/α)
    return 1/(2*α*pi) * ϵ^(1 + (1 - β)/α)*exp(ϵ^(1/α) * cos(ϕ/α)) * (cos(ω) + im*sin(ω))/(ϵ*exp(im*ϕ)-z)
end

Pint(α, β, ϵ, z) = ourquadgk(ϕ -> P(α, β, ϵ, ϕ, z), -α*pi, α*pi)
Pint(α, β, z) = Pint(α, β, 1, z)

function Kint(α, β, a, χ0, z)
    return ourquadgk(χ -> K(α, β, χ, z), a, χ0)
end

function K(α, β, r, z)
    den = (r^2 - 2*r*z*cos(α*pi) + z^2)
    return 1/(α*pi) * r^((1-β)/α) * exp(-r^(1/α))*( r*sin(pi*(1-β)) - z*sin(pi*(1 - β + α)))/den
end

Kint(α, β, χ0, z) = Kint(α, β, 0, χ0, z)

mpow(x::Complex, y) = x^y
mpow(x::Real, y) = x >= 0 ? x^y : Complex(x, 0)^y

# The following gymnastics are to get around the type-unstable mpow
function sum2(α, β, z::Real, k0)
    return z > 0 ? sum2_pos(α, β, z, k0) : sum2_neg(α, β, z, k0)
end
sum2(α, β, z::Complex, k0) = sum2_neg(α, β, z, k0)

for funcname in (:sum2_neg, :sum2_pos)
    local zarg
    local stype
    if funcname == :sum2_neg
        zarg = :( (z + Complex(0,0) ))
        stype = :( typeof(z + Complex(0,0) ) )
    else
        zarg = :z
        stype = :( typeof(z) )
    end
    @eval function ($funcname)(α, β, z, k0)
        s::($stype) = zero($stype)
        for k=1:k0
            arg = β - α * k
            if !( round(arg) == arg && arg < 0)
                s += ($zarg)^(-k) / gamma(arg)
            end
        end
        return s
    end
end

function choosesum(α, β, z, ρ)
    k0 = floor(Int, -log(ρ)/log(abs(z)))
    if abs(angle(z)) < pi*α/4 + 1//2 * min(pi, pi*α)
        return 1/α * z^((1-β)/α) * exp(z^(1/α)) - sum2(α, β, z, k0)
    else
        return -sum2(α, β, z, k0)
    end
end

function mittleffsum(α, β, z)
    k0::Int = floor(Int, α) + 1
    s = zero(z)
    for k=0:(k0 - 1)
        s += mittleff(α / k0, β, mpow(z, (1 // k0)) * exp(2pi * im * k / k0))
    end
    return s / k0
end

# mittleffsum sometimes passes complex z with magnitude nearly equal to 1
# Then k0 is astronomically high. Because the gamma function is evaluated by a call to gsl,
# the loop is not interruptible.
# So, we compute k0 using only the real part of k0.
# The effect on the accuracy is not controlled.
function mittleffsum2(α, β, z, ρ)
    zr = real(z)
    k0 = max(ceil(Int, (1-β)/α), ceil(Int, log(ρ*(1-abs(zr)))/log(abs(zr))))
    s = zero(z)
    for k=0:k0
        s += z^k/gamma(β+α*k)
    end
    return s
end


function mittleffints(α, β, z, ρ)
    az = abs(z)
    ab = abs(β)
    χ0 = β >= 0 ?
      max(1, 2*az, (-log(pi*ρ/6))^α) :
      max((ab+1)^α, 2*az, (-2*log( pi*ρ/(6*(ab+2)*(2*ab)^ab)))^α)
    aaz = abs(angle(z))
    if aaz > α * pi
        if β <= 1
            return Kint(α, β, χ0, z)
        else
            return Kint(α, β, 1, χ0, z) + Pint(α, β, z)
        end
    elseif aaz < pi*α
        if β <= 1
            return Kint(α, β, χ0, z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        else
            return Kint(α, β, az/2, χ0, z) + Pint(α, β, (az)/2, z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        end
    else
        return Kint(α, β, (az+1)/2, χ0, z) + Pint(α, β, (az+1)/2, z)
    end
end

"""
    mittlefferr(α, z, ρ)

Compute mittlefferr(α,1,z,ρ).
"""
mittlefferr(α, z, ρ) = mittlefferr(α, 1, z, ρ)


"""
    mittlefferr(α, β, z, ρ)

Compute the Mittag-Leffler function at `z` for parameters `α, β` with
accuracy `ρ`.
"""
function mittlefferr(α, β, z, ρ::Real)
    ρ > 0 || throw(DomainError(ρ))
    _mittlefferr(α, β, z, ρ)
end

_mittlefferr(α::Real, β::Real, z::Real, ρ::Real) = real(_mittleff(α, β, z, ρ))
_mittlefferr(α::Real, β::Real, z::Complex, ρ::Real) = _mittleff(α, β, z, ρ)

# The second definition would work for both complex and real
myeps(x) = x |> one |> float |> eps
myeps(x::Complex) =  x |> real |> myeps

"""
    mittleff(α, β, z)

Compute the Mittag-Leffler function at `z` for parameters `α, β`.
"""
mittleff(α, β, z) = _mittleff(α, β, float(z))
mittleff(α, β, z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, β, float(z))

#mittleff(α, β, z) = _mittlefferr(α, β, z, myeps(z))

"""
    mittleff(α, z)

Compute `mittleff(α, 1, z)`.
"""
mittleff(α, z) = _mittleff(α, 1, z)
#mittleff(α, z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, float(z))

mittleff(α, vec::Vector) = map(x -> mittleff(α, 1, x), vec)
mittleff(α, β, vec::Vector) = map(x -> mittleff(α, β, x), vec)
mittleff(α, β, γ, vec::Vector) = map(x -> mittleff(α, β, γ, x), vec)

function _mittleff_special_beta_one(α, z)
    z == 0 && return one(z)
    if α == 1/2 && abs(z) < 10
        return exp(z^2) * erfc(-z)
    end
    α == 0 && return 1/(1-z)  # FIXME needs domain check
    α == 1 && return exp(z)
    zc = isa(z, Real) && z < 0 ? complex(z) : z # Julia sometimes requires explicit complex numbers for efficiency
    α == 2 && return cosh(sqrt(zc))
    α == 3 && return (1//3) * (exp(zc^(1//3)) + 2 * exp(-zc^(1//3)/2) * cos(sqrt(convert(typeof(zc),3))/2 * zc^(1//3)))
    α == 4 && return (1//2) * (cosh(zc^(1//4)) + cos(zc^(1//4)))
    return nothing
end

function _mittleff_slow_with_eps(α, β, z, ρ)
    az = abs(z)
    az < 1 && return mittleffsum2(α, β, z, ρ)
    az > floor(10 + 5*α) && return choosesum(α, β, z, ρ)
    return mittleffints(α, β, z, ρ)
end

_mittleff(α::Real, β::Real, z::Real) = real(_mittleff0(α, β, z))
_mittleff(α, β, z) = _mittleff0(α, β, z)

function _mittleff0(α, β, z)
    if β == 1
        res = _mittleff_special_beta_one(α,z)
        res != nothing && return res
    end
    z == 0 && return 1/gamma(β)
    if β == 2
        α == 1 && return (exp(z) - 1)/z
        zc = isreal(z) && z < 0 ? complex(z) : z
        α == 2 && return sinh(sqrt(zc))/sqrt(zc)
    end
    1 < α && return mittleffsum(α,β,z)
    ρ = myeps(z)
    return _mittleff_slow_with_eps(α, β, z, ρ)
end

function _mittleff(α, β, z, ρ)
    if β == 1
        res = _mittleff_special_beta_one(α, z)
        res != nothing && return res
    end
    z == 0 && return 1/gamma(β)
    α == 1 && β == 1 && return(exp(z))
    1 < α && return mittleffsum(α, β, z)
    return _mittleff_slow_with_eps(α, β, z, ρ)
end

"""
    mittleffderiv(α, β, z)

Compute the derivative of the Mittag-Leffler function at `z` for parameters `α, β`.
"""
function mittleffderiv(α, β, z)
    # Derivative of Mittag Leffler function WRT to main argument Z.
    # Take q = 0.5. Paper requires |z| <= q < 1.
    q = 1//2
    #case 1, small z
    if abs(z) <= q
        ω = α + β - 3//2
        D = α^2 - 4*α*β + 6*α + 1
        #k1
        if α>1
            k₁ = ((2 - α - β)/(α - 1)) + 1
        elseif α > 0 && α <= 1 && D <= 0
            k₁ = ((3 - α - β)/α) + 1
        else
            k₁ = maximum([((3-α-β)/α) + 1, ((1-2*ω*α + sqrt(D))/(2*(α^2)))+1])
        end
        k₀ = maximum([k₁, log(myeps(z)*(1-abs(z)))/log(abs(z))])
        k₀ = ceil(Int, k₀) #take ceiling (not specified in paper whether floor or ceiling)
        out = zero(z)
        for k in 0:k₀
            out = out + ((k + 1)*z^k)/(gamma(α + β + α*k))
        end
    #case 2, larger z
    else
        out = (mittleff(α, β-1, z) - (β - 1)*mittleff(α, β, z))/(α*z)
    end
    return out
end

"""
    mittleffderiv(α, z)

Compute mittleffderiv(α, 1, z)
"""
mittleffderiv(α, z) = mittleffderiv(α, 1, z)

"""
    mittleff(α, β, γ, z)

Compute three-parametric mittleff(α, β, γ, z).

```tex
@article{2015,
   title={Numerical Evaluation of Two and Three Parameter Mittag-Leffler Functions},
   volume={53},
   ISSN={1095-7170},
   url={http://dx.doi.org/10.1137/140971191},
   DOI={10.1137/140971191},
   number={3},
   journal={SIAM Journal on Numerical Analysis},
   publisher={Society for Industrial & Applied Mathematics (SIAM)},
   author={Garrappa, Roberto},
   year={2015},
   month={Jan},
   pages={1350–1369}
}
```
"""
function mittleff(alpha, beta, gamma, z)
    log_epsilon = log(10^(-15))

    E=0
    if abs(z) < 1e-15
        E=1/gamma(beta)
    else
        E=LTIversion(1, z, alpha, beta, gamma, log_epsilon)
    end
    return E
end

function LTIversion(t, lambda, alpha, beta, gama, log_epsilon)
    theta = angle(lambda)
    kmin = ceil(Int, -alpha/2 - theta/2/pi)
    kmax = floor(Int, alpha/2 - theta/2/pi)
    k_vett = collect(kmin : kmax)
    s_star = abs(lambda)^(1/alpha) * exp.(im*(theta.+2*k_vett*pi)/alpha)

    phi_s_star = (real(s_star).+ abs.(s_star))/2
    index_s_star = sortperm(phi_s_star)
    phi_s_star = sort(phi_s_star)
    s_star = s_star[index_s_star]
    index_save = phi_s_star .> 1e-15

    s_star = s_star[index_save]
    phi_s_star = phi_s_star[index_save]
    s_star = [0; s_star]
    phi_s_star = [0; phi_s_star]
    J1 = length(s_star)
    J = J1 - 1
    p = [max(0,-2*(alpha*gama-beta+1)) ; ones(J)*gama]
    q = [ones(J)*gama ; +Inf]
    phi_s_star = [phi_s_star; +Inf]
    admissible_regions = findall((phi_s_star[1:end-1] .< (log_epsilon - log(eps()))/t) .& (phi_s_star[1:end-1] .< phi_s_star[2:end]))
    JJ1 = admissible_regions[end]
    mu_vett = ones(ComplexF64, JJ1)*Inf
    N_vett = ones(ComplexF64, JJ1)*Inf
    h_vett = ones(ComplexF64, JJ1)*Inf
    find_region = 0

    while find_region==0
        for j1 = admissible_regions
            muj=0; hj=0; Nj=0
            if j1 < J1
                j1=Int64(j1)
                (muj, hj, Nj) = OptimalParam_RB(t, phi_s_star[j1], phi_s_star[j1+1], p[j1], q[j1], log_epsilon) ;
            else
                (muj, hj, Nj) = OptimalParam_RU(t,phi_s_star[j1],p[j1],log_epsilon) ;
            end
            mu_vett[j1] = muj; h_vett[j1] = hj; N_vett[j1] = Nj
        end
        if minimum(real.(N_vett[:])) > 200
            log_epsilon = log_epsilon + log(10)
        else
            find_region = 1
        end
    end

    (N, iN) = findmin(real.(N_vett)); mu = mu_vett[iN]; h = h_vett[iN];

    k = collect(-N : N)
    u = h.*k
    z = mu*(1*im*u.+1).^2
    zd = -2*mu*u .+ 2*mu*1*im
    zexp = exp.(z*t)
    F = z.^(alpha*gama-beta)./(z.^alpha .- lambda).^gama.*zd
    S = zexp.*F
    Integral = h*sum(S)/2/pi/im

    ss_star = s_star[(iN[1]+1):end]
    Residues = sum(1/alpha*(ss_star).^(1-beta).*exp.(t*ss_star))

    E = Integral + Residues
    if isreal(lambda) 
        E = real(E)
    end
end

function OptimalParam_RU(t, phi_s_star_j, pj, log_epsilon)
    sq_phi_s_star_j = sqrt(phi_s_star_j)
    if phi_s_star_j > 0
        phibar_star_j = phi_s_star_j*1.01
    else
    phibar_star_j = 0.01
    end
    sq_phibar_star_j = sqrt(phibar_star_j)
    
    f_min=1
    f_max=10
    f_tar=5
    
    stop=0
    sq_muj=0
    A=0
    Nj=0
    while stop==0
        phi_t = phibar_star_j*t
        log_eps_phi_t = log_epsilon/phi_t
        Nj=Complex(ceil(real(phi_t/pi*(1-3*log_eps_phi_t/2+sqrt(Complex(1-2*log_eps_phi_t))))), ceil(imag(phi_t/pi*(1-3*log_eps_phi_t/2+sqrt(Complex(1-2*log_eps_phi_t))))))
        A=pi*Nj/phi_t
        sq_muj = sq_phibar_star_j*abs(4-A)/abs(7-sqrt(1+12*A))
        fbar=((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)^(-pj)
        stop=(pj<1e-14) || (f_min<fbar && fbar<f_max)
        if stop==0
            sq_phibar_star_j = f_tar^(-1/pj)*sq_muj+sq_phi_s_star_j
            phibar_star_j=sq_phibar_star_j^2
        end
    end
    muj=sq_muj^2
    hj=(-3*A-2+2*sqrt(1+12*A))/(4-A)/Nj
    
    log_eps=log(eps())
    threshold = (log_epsilon-log_eps)/t
    if muj>threshold
        if abs(pj)<1e-14
            Q=0
        else
            Q=f_tar^(-1/pj)*sqrt(muj)
        end
        phibar_star_j=(Q+sqrt(phi_s_star_j))^2
        if phibar_star_j<threshold
            w=sqrt(log_eps/(log_eps-log_epsilon))
            u=sqrt(-phibar_star_j*t/log_eps)
            muj = threshold
            Nj=ceil(w*log_epsilon/2/pi/(u*w-1))
            hj=sqrt(log_eps/(log_eps-log_epsilon))/Nj
        else
            Nj=+Inf
            hj=0
        end
    end
    return muj, hj, Nj
end

function OptimalParam_RB(t, phi_s_star_j, phi_s_star_j1, pj, qj, log_epsilon)
    log_eps = -36.043653389117154
    fac = 1.01
    conservative_error_analysis = 0
    
    f_max = exp(log_epsilon - log_eps)
    
    sq_phi_star_j = sqrt(phi_s_star_j)
    threshold = 2*sqrt((log_epsilon - log_eps)/t)
    sq_phi_star_j1 = min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j)
    
    if pj < 1.0e-14 && qj < 1.0e-14
        sq_phibar_star_j = sq_phi_star_j
        sq_phibar_star_j1 = sq_phi_star_j1
        adm_region = 1
    end
    
    if pj < 1.0e-14 && qj >= 1.0e-14
        sq_phibar_star_j = sq_phi_star_j
        if sq_phi_star_j > 0
            f_min = fac*(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j))^qj ;
        else
            f_min = fac
        end
        if f_min < f_max
            f_bar = f_min + f_min/f_max*(f_max-f_min)
            fq = f_bar^(-1/qj)
            sq_phibar_star_j1 = (2*sq_phi_star_j1-fq*sq_phi_star_j)/(2+fq) ;
            adm_region = 1
        else
            adm_region = 0 ;
        end
    end
    
    if pj >= 1.0e-14 && qj < 1.0e-14
        sq_phibar_star_j1 = sq_phi_star_j1
        f_min = fac*(sq_phi_star_j1/(sq_phi_star_j1-sq_phi_star_j))^pj
        if f_min < f_max
            f_bar = f_min + f_min/f_max*(f_max-f_min)
            fp = f_bar^(-1/pj)
            sq_phibar_star_j = (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp)
            adm_region = 1
        else
            adm_region = 0
        end
    end
    
    if pj >= 1.0e-14 && qj >= 1.0e-14
        f_min = fac*(sq_phi_star_j+sq_phi_star_j1)/(sq_phi_star_j1-sq_phi_star_j)^max(pj,qj) ;
        if f_min < f_max
            f_min = max(f_min,1.5)
            f_bar = f_min + f_min/f_max*(f_max-f_min)
            fp = f_bar^(-1/pj)
            fq = f_bar^(-1/qj)
            if conservative_error_analysis==0
                w = -phi_s_star_j1*t/log_epsilon
            else
                w = -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t)
            end
            den = 2+w - (1+w)*fp + fq
            sq_phibar_star_j = ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den ;
            sq_phibar_star_j1 = (-(1+w)*fq*sq_phi_star_j + (2+w-(1+w)*fp)*sq_phi_star_j1)/den ;
            adm_region = 1
        else
            adm_region = 0
        end
    end
    if adm_region==1
        log_epsilon = log_epsilon - log(f_bar)
        if conservative_error_analysis==0
            w = -sq_phibar_star_j1^2*t/log_epsilon
        else
            w = -2*sq_phibar_star_j1^2*t/(log_epsilon-sq_phibar_star_j1^2*t) ;
        end
        muj = (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2 ;
        hj = -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)/((1+w)*sq_phibar_star_j + sq_phibar_star_j1) ;
        Nj = Complex(ceil(real(sqrt(Complex(1-log_epsilon/t/muj))/hj)), ceil(imag(sqrt(Complex(1-log_epsilon/t/muj))/hj)))
    else
        muj = 0 ; hj = 0 ; Nj = +Inf
    end
    return muj, hj, Nj
end


function mlds(z,al,be,k)
    max_gamma_arg = 171.624
    Jmax = floor(Int, (max_gamma_arg - be)/al)
    G = gamma.(al.*(collect(0:Jmax)).+be)
    jj = collect(k:Jmax)
    f = ones(size(jj))
    for l = 1 : k
        f = f.*(jj.-l.+1)
    end
    c = f./G[k+1:Jmax+1]
    E = zeros(size(z)) ; Err_Round = zeros(size(z))
    Err_Round1 = copy(Err_Round)
    Err_Round2 = copy(Err_Round)
    for n = 1 : length(z)
        if abs(z[n]) < eps()
            E[n] = factorial(k)/G[k+1]
        else
            sum_arg = c.*z[n].^(jj.-k)
            abs_sum_arg = abs.(sum_arg)
            i_abs_sum_arg = abs_sum_arg .> eps()/2
            if any(i_abs_sum_arg) == 0
                i_abs_sum_arg[1] = 1
            end
            abs_sum_arg = abs_sum_arg[i_abs_sum_arg]
            sum_arg = sum_arg[i_abs_sum_arg]
            i_abs_sum_arg = sortperm(abs_sum_arg)
            abs_sum_arg = sort(abs_sum_arg)
            sum_arg = sum_arg[i_abs_sum_arg]
            if length(sum_arg) == 1
                E[n] = sum_arg
            else
                S = cumsum(sum_arg)
                S = S[2:end]
                E[n] = S[end]
            end
            J = length(sum_arg) - 1
            JJ = [J; collect(J:-1:1)]
            Err_Round1[n] = sum(JJ.*abs_sum_arg)*eps()
            if length(sum_arg) == 1
                Err_Round2[n] = Err_Round1[n]
            else
                Err_Round2[n] = sum(abs.(S))*eps()
            end
            Err_Round[n] = exp((log(Err_Round1[n])+log(Err_Round2[n]))/2)
        end
    end
    return E, Err_Round
end

function mldr(t,s,alpha,beta,k)
    omega = zeros(k+2)
    omega[1] = 1
    jj = collect(0:k)
    p = (alpha.-jj)
    pr = 1
    for j = 1:k+1
        pr = pr*p[j] ;
        omega[j+1] = pr/factorial(j) ;
    end
    Hk = zeros(k+1) ; Hk[1] = 1 ;
    for j = 1 : k
        ll = collect(1:j)
        Hk[j+1] = -1 ./alpha.*sum(omega[ll.+2].*(k.*ll./j.+1).*Hk[j.-ll.+1])  ;
    end
    ck = zeros(k+1)
    for j = 0:k
        temp = 0
        for l = 0 : k-j
            if l == 0
                p = 1 ;
            else
                ll = 0:l-1 ;
                p = prod(alpha.-beta.-ll) ;
            end
            temp = temp .+ p*Hk[k-j-l+1]./factorial(l) ;
        end
        ck[j+1] = temp/factorial(j) ;
    end
    tep = Polynomial(ck)
    result=tep.(s)
    R = 1 ./alpha^(k+1)*exp.(t.*s).*(s//1)^(1-alpha.*k.-beta).*result
    return R
end
#=
"""
```tex
@article{Garrappa2018ComputingTM,
  title={Computing the Matrix Mittag-Leffler Function with Applications to Fractional Calculus},
  author={Roberto Garrappa and Marina Popolizio},
  journal={Journal of Scientific Computing},
  year={2018},
  volume={77},
  pages={129-153}
}
```
"""
function mldlt(z,alpha,beta,k)
    log_epsilon = log(10^(-15)) ; 
    E = zeros(size(z)) ;  
    for j = 1 : length(z)
        if abs(z[j] < 1.0e-14)
            E[j] = factorial(k)/gamma(alpha*k+beta)
        else
            E[j] = LTIversion(1,z[j],alpha,beta,k,log_epsilon) ;
        end
        if isreal(z[j])
            E[j] = real(E[j])
        end 
    end
    return E
end

function mld(z, alpha, beta, k)
    Ek = zeros(size(z))
    tau = 1.0e-14 ; 
    max_gamma_arg = 171.624 ;
    Jmax = floor((max_gamma_arg - beta)/alpha) ;
    z_abs_max = (tau*gamma(alpha*Jmax+beta)/prod(Jmax.-collect(0:k-1))).^(1/(Jmax-k)) ;
    i_z_se = abs.(z).<z_abs_max ;
    z_se = z[i_z_se] ;
    E_se = zeros(size(z)) ; Err_Round = zeros(size(z)) ; e = ones(size(z)) ;
    (E_se[i_z_se], Err_Round[i_z_se]) = mlds(z_se,alpha,beta,k) ;
    i_z_se_accept = (Err_Round./(e.+abs.(E_se)) .< tau) .& ~(E_se==0)  ;
    Ek[i_z_se_accept] = E_se[i_z_se_accept]
    i_ze_lt = .!i_z_se_accept ;
    if k <= 3
        p = 0 
    elseif k <= 7
        p = 1
    else
        p = 2
    end 
    c = zeros(k-p+1, k-p+1) ;
    c[1, 1] = 1 ;
    for kk = 1 : k-p
        c[kk+1,1] = (1-alpha*(kk-1)-beta)*c[kk,1]
        for j = 1 : kk-1
            c[kk+1,j+1] = c[kk,j] + (j+1-alpha*(kk-1)-beta)*c[kk,j+1]
        end
        c[kk+1,kk+1] = 1 ;
    end
    for j = 0 : k-p
        if abs(c[k-p+1,j+1]) > 1.0e-14
            Ek[i_ze_lt] = Ek[i_ze_lt] .+ c[k-p+1,j+1]*mldlt(z[i_ze_lt], alpha, (k-p)*alpha+beta-j, p) ;
        end
    end
    Ek[i_ze_lt] = Ek[i_ze_lt]./alpha^(k-p) ;
    return Ek
end

function mittleff(A, alpha, beta)
    fun(z, k) = mld(z, alpha, beta, k)
    E=funm(A)
end
#FIXME: Didn't see any matrix function support in Julia😟.
# Yingbo has a implementation: https://github.com/YingboMa/Funm.jl
=#