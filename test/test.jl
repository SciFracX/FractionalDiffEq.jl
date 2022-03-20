using FractionalDiffEq, SpecialFunctions
using Test

@testset "Test Diethelm PECE algorithms" begin
    fun(x, y) = 1-y
    prob = SingleTermFODEProblem(fun, 1.8, 0, 5)
    result = solve(prob, 0.01, PECE())
    tspan = collect(0:0.01:5)

    target = []

    for i in 0:0.01:5
        push!(target, i^1.8*mittleff(1.8, 2.8,-i^1.8))
    end

    @test isapprox(result, target; atol=1)
end

@testset "Test Matrix discrete method" begin
    fun(x, y) = 1-y
    tspan=collect(0.01:0.01:5)
    target = []
    for i in 0.01:0.01:5
        push!(target, i^0.5*mittleff(0.5, 1.5,-i^0.5))
    end

    singleprob = MultiTermsFODEProblem([1, 1], [0.5, 0], 1, 5)
    result = solve(singleprob, 0.01, FODEMatrixDiscrete())

    @test isapprox(result, target; atol=1)

    ########### Yet another test case ##########
    yafun(x, y) = 1-y

    #Analytical solution
    yatarget = []

    #MittagLeffler.jl doesn't support array argument
    for i in 0.01:0.01:20
        push!(yatarget, i^1.8*mittleff(1.8,2.8,-i^1.8))
    end

    highsingleprob = MultiTermsFODEProblem([1, 1], [1.8, 0], 1, 20)
    yaresult = solve(highsingleprob, 0.01, FODEMatrixDiscrete())

    @test isapprox(yaresult, yatarget; atol=1)

end

@testset "Test Matrix discrete method for FPDE" begin
    result = solve(0.7, 1.8, 1, 1, 3, 2, FPDEMatrixDiscrete())
    @test result ≈ [0.0  0.0; 0.0  0.0644409; 0.0  0.0]
end

@testset "Test FiniteDiffEx method" begin
    K = 1
    fdorder = 1.9
    dx = pi/2
    dt = 0.5
    n = 2
    xStart = 0
    xEnd = pi
    u0t = 0
    uendt = 0
    u0(x) = sin(x)

    U=solve(fdorder, dx, dt, xStart, xEnd, n, K, u0t, uendt, u0, FiniteDiffEx())
    @test isapprox(U, [ 0.0  1.0  1.22465e-16; 0.0   0.793379    0.0; 0.0   0.43766     0.0; 0.0   0.00221205  0.0; 0.0 -0.42797 0.0]; atol=1e-3)
end

@testset "Test FiniteDiffIm method" begin
    K = 1
    α = 0.5
    dx = pi/2
    dt = 0.5
    n = 2
    xStart = 0
    xEnd = pi
    u0t = 0
    uendt = 0
    u0(x) = sin(x)

    U=solve(α, dx, dt, xStart, xEnd, n, K, u0t, uendt, u0, FiniteDiffIm())
    @test isapprox(U, [ 0.0  1.0       1.22465e-16
    0.0  0.663152  0.0
    0.0  0.532299  0.0
    0.0  0.459938  0.0
    0.0  0.412321  0.0]; atol=1e-3)
end

@testset "Test Caputo Discrete for Advection Diffusion equation" begin
    fx0(x) = 0
    function fgz(q)
        x = q[1, 2];t=q[1, 1];a=q[1, 3]
        f = exp(x)*t^(6-a)/gamma(7-a)*720
        return f
    end
    function f0t(k)
        x=k[1, 1]
        a=k[1, 2]
        f=x^6
        return f
    end
    function flt(k)
        x = k[1,1];a=k[1,2]
        f = exp(1)*x^6
        return f
    end
    result = solve(3, 0.5, 1, 2, 1, 20, 1, 1, fx0, fgz, f0t, flt, ADV_DIF())
    @test isapprox(result, [0.0  0.0171843  1.06757
    0.0  0.0187449  1.1369
    0.0  0.0203056  1.20805
    0.0  0.0218648  1.2811
    0.0  0.0234204  1.35608
    0.0  0.0249696  1.43307
    0.0  0.0265091  1.51209
    0.0  0.028035   1.59319
    0.0  0.0295425  1.67637
    0.0  0.0310258  1.76166
    0.0  0.0324786  1.84904
    0.0  0.0338932  1.93849
    0.0  0.0352607  2.02997
    0.0  0.0365713  2.12342
    0.0  0.0378133  2.21874
    0.0  0.0389738  2.31582
    0.0  0.0400378  2.41452
    0.0  0.0409888  2.51465
    0.0  0.0418076  2.61599]; atol=1e-3)
end

@testset "Test Closed Form method" begin
    t=collect(0:0.002:10);
    rightfun(x) = sin(x)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3], t)
    result=solve(prob, ClosedForm())
end

@testset "Test ClosedFormHankelM method" begin
    t = collect(0:0.5:1);
    rightfun(x)=sin(x^2)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3], t)
    result = solve(prob, ClosedFormHankelM())

    @test result≈[0.0; 0.08402140107687359; 0.3754974742112727]
end

@testset "Test GLWithMemory method" begin
    h=0.5
    alpha = [0.99, 0.99, 0.99]
    x0 = [1, 0, 1]
    tf=1
    function f(t, x, y, z, k)
        a, b, c = 10, 28, 8/3
        if k == 1
            return a*(y-x)
        elseif k == 2
            return x*(b-z)-y
        elseif k == 3
            return x*y-c*z
        end
    end
    prob = FODESystem(f, alpha, x0)
    result = solve(prob, h, tf, GLWithMemory())
    @test isapprox(result, [1.0 0.0 1.0
    -4.04478 13.5939 -0.352607       
    84.8074 -51.1251 -27.5541]; atol=1e-4)
end

@testset "Test GL method" begin
    fun(x, y) = 1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 1)
    result = solve(prob, 0.1, GL())

    @test isapprox(result, [0.0
    0.31622776601683794
    0.3743416490252569
    0.42454983788325495
    0.4608051796660425
    0.4897286932245971
    0.5136308879844076
    0.5339402943673064
    0.5515444532932976
    0.5670394859746068
    0.5808448788127054]; atol=1e-3)
end

@testset "Test DelayPECE method" begin
    function ϕ(x)
        if x == 0
            return 19.00001
        else
            return 19.0
        end
    end
    
    function f(t, y, ϕ)
        return 3.5*y*(1-ϕ/19)
    end
    
    h = 0.5
    α = 0.97
    τ = 0.8
    T = 1
    fddeprob = FDDEProblem(f, ϕ, α, τ, T)
    V, y = solve(fddeprob, h, DelayPECE())

    @test V≈[19.0, 19.0, 1.0]
    @test y≈[19.00001, 19.00001, 37.274176448220274]
end

@testset "Test DelayPI method" begin
    function ϕ(x)
        if x == 0
            return 19.00001
        else
            return 19.0
        end
    end
    function f(t, y, ϕ)
        return 3.5*y*(1-ϕ/19)
    end
    prob = FDDEProblem(f, ϕ, 0.97, 0.8, 2, 0)
    result = solve(prob, 0.5, DelayPI())
    @test result≈[19.00001, 19.00001, 19.00001, 18.99999190949352, 18.99997456359874]
end

@testset "Test DelayABM method" begin
    h=0.5; T=5; α=0.97; τ=2; q=0.5
    f(t, ϕ, y) = 2*ϕ/(1+ϕ^9.65)-y
    prob = FDDEProblem(f, q, α, τ, T)
    x, y=solve(prob, h, DelayABM())

    @test isapprox(x, [1.078559863692747, 1.175963999045738, 1.1661317460354588, 1.128481756921719, 1.0016061526083417, 0.7724564325042358, 0.5974978685646778]; atol=1e-3)
    @test isapprox(y, [0.8889787467894421, 0.9404487875504524, 0.9667449499617093, 0.9803311436135411, 1.078559863692747, 1.175963999045738, 1.1661317460354588]; atol=1e-3)
end

@testset "Test Nonlinear method" begin
    using FractionalDiffEq

    function chua(t, x, k)
        a = 10.725
        b = 10.593
        c = 0.268
        m0 = -1.1726
        m1 = -0.7872

        if k == 1
            f = m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))
            y = a*(x[2]-x[1]-f)
            return y
        elseif k == 2
            y = x[1]-x[2]+x[3]
            return y
        elseif k == 3
            y = -b*x[2]-c*x[3]
            return y
        end
    end

    α = [0.93, 0.99, 0.92];
    x0 = [0.2; -0.1; 0.1];
    h = 0.1;
    prob = FODESystem(chua, α, x0)
    tn = 0.5;
    result = solve(prob, h, tn, NonLinearAlg())

    @test isapprox(result, [0.2 -0.1 0.1
    0.11749    -0.0675115   0.182758       
    0.063749   -0.0357031   0.215719       
    0.0394769  -0.00641779  0.210729       
    0.0458196   0.019928    0.175057       
    0.0843858   0.0438356   0.113762]; atol=1e-3)
end

@testset "Test FDDE Matrix Form method" begin
    limit=100
    t0=0
    T=1
    tau=3.1416
    h=0.5
    alpha=0.4
    function testx0(t)
        return [sin(t)*cos(t); sin(t)*cos(t); cos(t)^2-sin(t)^2; cos(t)^2-sin(t)^2]
    end
    A=[0 0 1 0; 0 0 0 1; 0 -2 0 0; -2 0 0 0]
    B=[0 0 0 0; 0 0 0 0; -2 0 0 0; 0 -2 0 0]
    f=[0; 0; 0; 0]
    
    result=solve(limit, alpha, A, B, f, t0, testx0, T, tau, h, MatrixForm())

    @test isapprox(result, [0.139708   0.139708   0.96017    0.96017
    0.479462   0.479462   0.283662   0.283662
    0.378401   0.378401  -0.653644  -0.653644
   -0.07056   -0.07056   -0.989992  -0.989992
   -0.454649  -0.454649  -0.416147  -0.416147
   -0.420735  -0.420735   0.540302   0.540302
    0.0        0.0        1.0        1.0
    0.242142   0.242142   0.932342   0.932342
    0.401838   0.401838   0.60134    0.60134]; atol=1e-4)
end

@testset "Test Distributed Order Matrix method" begin    
    h = 0.5; t = collect(h:h:1);
    fun(t)=-0.1
    prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], [0, 1], t, fun)
    result = solve(prob, h, DOMatrixDiscrete())
    @test isapprox(result, [0; -0.06531153193179398]; atol=1e-3)
end

@testset "Test Fractional Differences Equations PECE method" begin
    differencefun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    prob = FractionalDifferenceProblem(differencefun, α, x0)
    t, y=solve(prob, T, h, PECEDifference())
    @test isapprox(y, [1.0, 2.5, 4.25, 5.125, 5.5625, 5.78125, 5.890625, 5.9453125, 5.97265625, 5.986328125, 5.9931640625]; atol=1e-3)
end

@testset "Test Fractional Integral Equations SpectralUltraspherical method" begin
    prob = FIEProblem([1, 1], [1, 0.5], 1)
    xx = LinRange(-1, 1, 5)
    sol=solve(prob, 20, xx, SpectralUltraspherical())

    @test isapprox(sol, [ 1.0
    0.5231565837302468
    0.42758357615580733
    0.37316567427801584
    0.3362040024463395]; atol=1e-3)
end

@testset "Test Product Integral Explicit method" begin
    fun(t, y)=1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 5)
    sol=solve(prob, 0.5, PIEx())

    @test isapprox(sol, [ 0.0
    0.5840920370824765
    0.6497603208539517
    0.6965874966593231
    0.7309421360673155
    0.7569810570702951
    0.7773192302138408
    0.7936279277715693
    0.8070043316131617
    0.8181895332990293]; atol=1e-4)
end

@testset "Test Product Integral Implicit method" begin
    fun(t, y)=1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 1)
    sol=solve(prob, 0.1, PIIm())

    @test isapprox(sol, [ 0.0
    0.504626504404032
    0.543454077521492
    0.5760954086164125
    0.6028545521580736
    0.6251322531561498
    0.6440148324633577
    0.6602762702077741
    0.6744705915156369
    0.6870025519994545]; atol=1e-4)
end

@testset "Test Product Integral Implicit Trapezoidal method" begin
    fun(t, y)=1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 1)
    sol=solve(prob, 0.1, PITrap())

    @test isapprox(sol, [0.0
    0.504626504404032
    0.5185925289496215
    0.5467126570149609
    0.571261200170596
    0.5923262209457972
    0.6105232343055992
    0.6264159719331727
    0.6404461835166287
    0.652951997319657]; atol=1e-4)
end