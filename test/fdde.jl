@testset "Test DelayPECE method with single constant lag" begin
    function h(p, t)
        if t == 0
            return 19.00001
        else
            return 19.0
        end
    end
    
    f(y, ϕ, p, t) = 3.5*y*(1-ϕ/19)
    
    dt = 0.2
    α = 0.97
    τ = [0.8]
    u0 = 19.00001
    tspan = (0.0, 1.0)
    fddeprob = FDDEProblem(f, α, u0, h, constant_lags = τ, tspan)
    sol = solve(fddeprob, DelayPECE(), dt=dt)

    @test isapprox(test_sol(sol)', [19.00001
    19.00001
    19.00001
    19.00001
    19.000006224428567
    18.999998984089164]; atol=1e-7)
end

@testset "Test DelayPECE method with single constant lag with variable order" begin
    function h(p, t)
        if t == 0
            return 19.00001
        else
            return 19.0
        end
    end
    
    f(y, ϕ, p, t) = 3.5*y*(1-ϕ/19)
    dt = 0.01
    alpha(t) = 0.99-(0.01/100)*t
    τ = [0.8]
    tspan = (0.0, 1.0)
    u0 = 19.00001
    fddeprob = FDDEProblem(f, [alpha], u0, h, constant_lags = τ, tspan)
    sol = solve(fddeprob, DelayPECE(), dt=dt)

    @test isapprox(test_sol(sol)'[end-10:end], [19.000006225431903
    19.00000586970618
    19.000005514290276
    19.000005159159087
    19.000004804291258
    19.000004449668406
    19.00000409527454
    19.000003741095618
    19.00000338711922
    19.000003033334277
    19.000002679730862]; atol=1e-7)
end

@testset "Test DelayPECE method with single time varying lag" begin
    function h(p, t)
        if t == 0
            return 19.00001
        else
            return 19.0
        end
    end
    
    function f(y, ϕ, p, t)
        return 3.5*y*(1-ϕ/19)
    end
    
    dt = 0.01
    alpha = 0.97
    u0 = 19.00001
    τ(t) = 0.8
    tspan = (0.0, 1.0)
    fddeprob = FDDEProblem(f, alpha, u0, h, constant_lags = [τ], tspan)
    sol = solve(fddeprob, DelayPECE(), dt=dt)

    @test isapprox(test_sol(sol)'[end-10:end], [19.000006018940965
    19.000005651668022
    19.000005285353794
    19.000004919919107
    19.000004555296744
    19.000004191428914
    19.00000382826542
    19.000003465762248
    19.000003103880502
    19.00000274258556
    19.00000238184641]; atol=1e-8)
end

@testset "Test Product Integral method" begin
    function ϕ(p, t)
        if t == 0
            return 19.00001
        else
            return 19.0
        end
    end
    function f(y, ϕ, p, t)
        return 3.5*y*(1-ϕ/19)
    end
    τ = [0.8]
    order = 0.97
    u0 = 19.00001
    tspan = (0.0, 2.0)
    dt = 0.5
    prob = FDDEProblem(f, order, u0, ϕ, constant_lags = τ, tspan)
    sol = solve(prob, DelayPIEX(), dt=dt)
    @test isapprox(test_sol(sol)', [19.00001, 19.00001, 19.00001, 18.99999190949352, 18.99997456359874]; atol=1e-3)
end
#DelayABM todo
#=
@testset "Test DelayABM method" begin
    dt=0.5
    tspan = (0.0, 5.0)
    α=0.97
    τ=[2.0]
    u0 = 0.5
    h(p, t) = 0.5
    delayabmfun(y, ϕ, p, t) = 2*ϕ/(1+ϕ^9.65)-y
    prob = FDDEProblem(delayabmfun, α, u0, h, constant_lags=τ, tspan)
    x, y=solve(prob, DelayABM(), dt=dt)

    @test isapprox(x, [1.078559863692747, 1.175963999045738, 1.1661317460354588, 1.128481756921719, 1.0016061526083417, 0.7724564325042358, 0.5974978685646778]; atol=1e-3)
    @test isapprox(y, [0.8889787467894421, 0.9404487875504524, 0.9667449499617093, 0.9803311436135411, 1.078559863692747, 1.175963999045738, 1.1661317460354588]; atol=1e-3)
end

@testset "Test DelayABM for FDDESystem" begin
    function EnzymeKinetics!(dy, y, ϕ, p, t)
        dy[1] = 10.5-y[1]/(1+0.0005*ϕ[4]^3)
        dy[2] = y[1]/(1+0.0005*ϕ[4]^3)-y[2]
        dy[3] = y[2]-y[3]
        dy[4] = y[3]-0.5*y[4]
    end
    function test_phi(p, t)
        return [60.0, 10.0, 10.0, 20.0]
    end
    tspan = (0.0, 6.0)
    q = [60.0, 10.0, 10.0, 20.0]; α = [0.95, 0.95, 0.95, 0.95]
    prob = FDDEProblem(EnzymeKinetics!, α, q, test_phi, constant_lags=[4], tspan)
    sol = solve(prob, DelayABM(), dt=0.1)
#=
    @test isapprox(sold, [ 56.3     11.3272   11.4284  22.4025
    56.3229  11.2293   11.4106  22.4194
    56.3468  11.141    11.3849  22.4332
    56.3731  11.0592   11.3533  22.4434
    56.4026  10.9824   11.3171  22.4495
    56.4358  10.9096   11.2773  22.4516
    56.4729  10.8398   11.2346  22.4493
    56.5144  10.7725   11.1896  22.4429
    56.5604  10.707    11.1429  22.4322
    56.6111  10.6431   11.0946  22.4175
    56.6668  10.5803   11.0452  22.3987
    56.7275  10.5186   10.9947  22.376
    56.7932  10.4577   10.9435  22.3495
    56.864   10.3975   10.8916  22.3194
    56.9398  10.3382   10.8392  22.2857
    57.0205  10.2795   10.7863  22.2488
    57.1059  10.2217   10.7332  22.2086
    57.1961  10.1647   10.6799  22.1653
    57.2908  10.1087   10.6265  22.1191
    57.3897  10.0537   10.5731  22.0702
    57.4928   9.99987  10.5198  22.0187]; atol=1e-3)
=#
    @test isapprox(sol.u, [60.5329  10.8225  10.6355  20.6245
    60.3612  10.9555  10.6638  20.6616
    60.1978  11.0679  10.7015  20.6988
    60.0407  11.1632  10.7456  20.7379
    59.8889  11.2442  10.7935  20.7792
    59.7419  11.3129  10.8437  20.8231
    59.5991  11.3711  10.8946  20.8697
    59.4603  11.4202  10.9452  20.9187
    59.3251  11.4613  10.9947  20.9701
    59.1935  11.4954  11.0425  21.0235
    59.0652  11.5234  11.0883  21.0786
    58.9401  11.5461  11.1315  21.1351
    58.818   11.5642  11.1722  21.1927
    58.6988  11.5781  11.21    21.2511
    58.5824  11.5883  11.245   21.3098
    58.4686  11.5954  11.2772  21.3687
    58.3575  11.5997  11.3065  21.4274
    58.2489  11.6016  11.3331  21.4857
    58.1427  11.6013  11.357   21.5433
    58.0389  11.5992  11.3783  21.6001
    57.9373  11.5953  11.3971  21.6557]'; atol=1e-3)
end

@testset "Test FDDE Matrix Form method" begin
    τ=3.1416; h=0.5; α=0.4; tspan = (0, 1)
    fx0(t) = [sin(t)*cos(t); sin(t)*cos(t); cos(t)^2-sin(t)^2; cos(t)^2-sin(t)^2]
    A=[0 0 1 0; 0 0 0 1; 0 -2 0 0; -2 0 0 0]
    B=[0 0 0 0; 0 0 0 0; -2 0 0 0; 0 -2 0 0]
    f=[0; 0; 0; 0]
    prob = FDDEMatrixProblem(α, τ, A, B, f, fx0, tspan)
    result=solve(prob, h, MatrixForm())

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

#=
@testset "Test DelayPECE method with multiple constant lags" begin
    α = 0.95; ϕ(x) = 0.5
    τ = [2, 2.6]
    fun(t, y, ϕ1, ϕ2) = 2*ϕ1/(1+ϕ2^9.65)-y
    prob = FDDEProblem(fun, ϕ, α, τ, 5)
    delayed, y = solve(prob, 0.5, DelayPECE())
    
    @test isapprox(y, [  0.5
    0.6920978716960207
    1.1225901176125745
    0.6230201580384503
    1.2802458016051212
    0.8184174623360388
    2.5359563307529904
   -1.6541509069477893
    6.108941265937457
   -5.480477100398328
   14.358155431964132]; atol=1e-3)

   @test isapprox(delayed, [0.5  0.5  0.5  0.5  0.5  1.12259   0.62302  1.28025   0.818417  2.53596  -1.65415
   0.5  0.5  0.5  0.5  0.5  0.605999  1.2225   0.491575  1.37261   0.47491   3.37398]; atol=1e-3)
end
=#
=#