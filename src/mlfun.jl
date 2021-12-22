## This code implements algorithms in
## Rudolfo Gorenflo, Joulia Loutchko and Yuri Loutchko,
## *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal, **(2002)**
"""
@ARTICLE{Gorenflo_computationof,
    author = {Rudolf Gorenflo and Joulia Loutchko and Yuri Luchko},
    title = {Computation of the Mittag-Leffler function Eα,β(z) and its derivative},
    journal = {FRACT. CALC. APPL. ANAL},
    year = {},
    pages = {2002}
}
"""
# This MittagLeffler function is modified from [John Lapeyre](https://github.com/jlapeyre)'s [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)
# Since Mittag Leffler functions is widely used in Fractional Differential Equation, so we decided to has the Mittag Leffler function build in.

import QuadGK: quadgk

function ourquadgk(f, a, b)
    (res, err) = quadgk(f, a, b; order=17)
    return res
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
        return - sum2(α, β, z, k0)
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
    mittlefferr(α,z,ρ)

Compute mittlefferr(α,1,z,ρ).
"""
mittlefferr(α, z, ρ) = mittlefferr(α, 1, z, ρ)


"""
    mittlefferr(α,β,z,ρ)

Compute the Mittag-Leffler function at `z` for parameters `α,β` with
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

Compute the Mittag-Leffler function at `z` for parameters `α,β`.
"""
mittleff(α, β, z) = _mittleff(α, β, float(z))
mittleff(α, β, z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, β, float(z))

#mittleff(α, β, z) = _mittlefferr(α,β,z,myeps(z))

"""
    mittleff(α, z)

Compute `mittleff(α, 1, z)`.
"""
mittleff(α, z) = _mittleff(α,1,z)
#mittleff(α, z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, float(z))

function _mittleff_special_beta_one(α,z)
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
    mittleffderiv(α,β,z)

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
    mittleffderiv(α,z)

Compute mittleffderiv(α,1,z)
"""
mittleffderiv(α, z) = mittleffderiv(α, 1, z)