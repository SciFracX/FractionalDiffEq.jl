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
    u0f(x) = sin(x)

    U=solve(fdorder, dx, dt, xStart, xEnd, n, K, u0t, uendt, u0f, FiniteDiffEx())
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
    u0f(x) = sin(x)

    U=solve(α, dx, dt, xStart, xEnd, n, K, u0t, uendt, u0f, FiniteDiffIm())
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