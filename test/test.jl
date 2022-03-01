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

    singleprob = MultiTermsFODEProblem([1, 1], [0.5, 0], 1)
    result = solve(singleprob, 0.01, 5, FODEMatrixDiscrete())

    @test isapprox(result, target; atol=1)

    ########### Yet another test case ##########
    yafun(x, y) = 1-y

    #Analytical solution
    yatarget = []

    #MittagLeffler.jl doesn't support array argument
    for i in 0.01:0.01:20
        push!(yatarget, i^1.8*mittleff(1.8,2.8,-i^1.8))
    end

    highsingleprob = MultiTermsFODEProblem([1, 1], [1.8, 0], 1)
    yaresult = solve(highsingleprob, 0.01, 20, FODEMatrixDiscrete())

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

    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3])

    result=solve(prob, t, ClosedForm())

end

@testset "Test ClosedFormHankelM method" begin
    t = collect(0:0.5:1);
    rightfun(x)=sin(x^2)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3])
    result = solve(prob, t, ClosedFormHankelM())

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
    x, y, z = solve(prob, h, tf, GLWithMemory())
    @test x≈[1.0, -4.044777750283594, 84.80744193501619]
    @test y≈[0.0, 13.593899925765704, -51.12509457411144]
    @test z≈[1.0, -0.3526074000756252, -27.554093040332816]
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
    fddeprob = FDDEProblem(f, ϕ, α, τ)
    V, y = solve(fddeprob, T, h, DelayPECE())

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
    prob = FDDEProblem(f, ϕ, 0.97, 0.8, 0)
    result = solve(prob, 2, 0.5, DelayPI())
    @test result≈[19.00001, 19.00001, 19.00001, 18.99999190949352, 18.99997456359874]
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
    B=[0 0 0 0; 0 0 0 0 ;-2 0 0 0; 0 -2 0 0]
    f=[0; 0; 0; 0]
    
    result=solve(limit, t0, T, tau, h, alpha, testx0, A, B, f, MatrixForm())

    @test isapprox(result, [ 0.139708  0.479462   0.378401  -0.07056   -0.454649  -0.420735  0.0  0.242142  0.401838
    0.139708  0.479462   0.378401  -0.07056   -0.454649  -0.420735  0.0  0.242142  0.401838
    0.96017   0.283662  -0.653644  -0.989992  -0.416147   0.540302  1.0  0.932342  0.60134
    0.96017   0.283662  -0.653644  -0.989992  -0.416147   0.540302  1.0  0.932342  0.60134]; atol=1e-4)
end

@testset "Test Distributed Order Matrix method" begin    
    h = 0.5; t = collect(h:h:1);
    fun(t)=-0.1
    prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], [0, 1], t, h, fun)
    result = solve(prob, DOMatrixDiscrete())
    @test isapprox(result, [0; -0.06531153193179398]; atol=1e-3)
end

@testset "Test Fractional Differences Equations PECE method" begin
    fun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    t, y=solve(fun, α, x0, T, h, PECEDifference())
    @test isapprox(y, [1.0, 2.5, 4.25, 5.125, 5.5625, 5.78125, 5.890625, 5.9453125, 5.97265625, 5.986328125, 5.9931640625]; atol=1e-3)
end