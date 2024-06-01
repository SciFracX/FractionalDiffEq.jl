@testset "Test AtanganaSeda method for FFODEProblem" begin
    α = 1
    β = 1
    h = 0.1
    tfinal = 0.5
    u0 = [-2, 1, -1]
    a = 10
    b = 28
    c = 8 / 3
    function fun(du, u, p, t)
        a = 10
        b = 28
        c = 8 / 3
        du[1] = a * (u[2] - u[1])
        du[2] = (b - u[3]) * u[1] - u[2]
        du[3] = u[1] * u[2] - c * u[3]
    end
    prob = FFMODESystem(fun, [α, β], u0, (0, tfinal))
    sol = solve(prob, h, AtanganaSeda())

    @test isapprox(sol.u,
        [-2.0 1.0 -9.35 31.0271 -160.86 637.491
         1.0 -4.9 3.125 -59.1272 190.243 -11690.4
         -1.0 -0.933333 -1.32833 -5.57208 -351.022 -5795.54];
        atol = 1e-2)
end

@testset "Test AtanganaSeda method for Caputo-Fabrizio sense variable order FFODE" begin
    alpha = 0.96
    h = 0.1
    tfinal = 2
    bet(t) = 0.01 + 0.01 * t
    u0 = [-0.2; 0.5; 0.2]
    function fun(du, u, p, t)
        gama = 10.814
        lambda = 14
        a = 0.2
        b = 0.15
        du[1] = gama * (u[2] - a * sin(2 * pi * b * u[1]))
        du[2] = u[1] - u[2] + u[3]
        du[3] = -lambda * u[2]
    end
    prob = FFMODESystem(fun, [alpha, bet], u0, (1, tfinal))
    result = solve(prob, h, AtanganaSeda())
    @test isapprox(result.u,
        [-0.2 0.381227 0.706487 0.705708 0.709295 0.712962 0.716751 0.720653 0.724658 0.72876 0.73295
         0.5 0.45 0.389684 0.388033 0.387055 0.386023 0.384947 0.383829 0.382669 0.38147 0.380232
         0.2 -0.5 -1.095 -1.09861 -1.10538 -1.11243 -1.11973 -1.12727 -1.13504 -1.14303 -1.15123];
        atol = 1e-3)
end

@testset "Test AtanganaSeda method for single term FFMODEProblem" begin
    fun(t, y) = -3 * t + 17
    u0 = -0.1
    h = 0.01
    α = 0.75
    β = 0.5
    tfinal = 0.1
    tspan = (0, tfinal)
    prob = FFMODESystem(fun, [α, β], u0, tspan)
    sol = solve(prob, h, AtanganaSeda())

    isapprox(sol.u,
        [-0.1, 0.06969999999999998, 0.23895000000000002, 1.5769475089199805,
            3.261575723419437, 4.805960383747965, 6.0083131151268985,
            6.860609872742969, 7.439542028856125, 7.834743031138744];
        atol = 1e-4)
end
