#=
@testset "Test Matrix discrete method" begin
    fun(x, y) = 1-y
    u0 = [0, 0]
    tspan=collect(0.01:0.01:5)
    target = []
    for i in 0.01:0.01:5
        push!(target, i^0.5*mittleff(0.5, 1.5,-i^0.5))
    end

    singleprob = MultiTermsFODEProblem([1, 1], [0.5, 0], 1, u0, (0, 5))
    result = solve(singleprob, 0.01, FODEMatrixDiscrete())

    @test isapprox(result.u, target; atol=1)

    ########### Yet another test case ##########
    yafun(x, y) = 1-y

    #Analytical solution
    yatarget = []

    #MittagLeffler.jl doesn't support array argument
    for i in 0.01:0.01:20
        push!(yatarget, i^1.8*mittleff(1.8,2.8,-i^1.8))
    end

    highsingleprob = MultiTermsFODEProblem([1, 1], [1.8, 0], 1, u0, (0, 20))
    yaresult = solve(highsingleprob, 0.01, FODEMatrixDiscrete())

    @test isapprox(yaresult.u, yatarget; atol=1)

end

@testset "Test Closed Form method" begin
    rightfun(x)=sin(x^2)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90 3], [1 0.3 0], [0, 0], (0, 5))
    sol = solve(prob, 0.5, ClosedForm())

    @test isapprox(sol.u, [ 0.0
    0.0854971609134175
    0.38289836071474537
    0.6454389110668579
    0.2783045139283708
    0.043460434414588106
    0.12348574062908478
    0.04623359232392658
   -0.07835912500702993
    0.26027660764161137
    0.2855329523586738]; atol=1e-3)
end

@testset "Test ClosedFormHankelM method" begin
    u0 = [0, 0, 0, 0]
    rightfun(x)=sin(x^2)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3], u0, (0, 1))
    result = solve(prob, 0.5, ClosedFormHankelM())

    @test isapprox(result.u, [0.0; 0.08402140107687359; 0.3754974742112727]; atol=1e-3)
end
=#
@testset "Test GL method for FODEProblem" begin
    alpha = [0.99, 0.99, 0.99]
    u0 = [1.0, 0.0, 1.0]
    tspan = (0.0, 1.0)
    function testf!(du, u, p, t)
        a, b, c = 10, 28, 8/3
        du[1] = a*(u[2]-u[1])
        du[2] = u[1]*(b-u[3])-u[2]
        du[3] = u[1]*u[2]-c*u[3]
    end
    prob = FODEProblem(testf!, alpha, u0, tspan)
    sol = solve(prob, GL(), dt=0.5)
    @test isapprox(test_sol(sol), [1.0  -4.04478    84.8074
    0.0  13.5939    -51.1251
    1.0  -0.352607  -27.5541]; atol=1e-4)
end

@testset "Test Nonlinear method" begin
    function chua!(du, x, p, t)
        a = 10.725; b = 10.593
        c = 0.268
        m0 = -1.1726
        m1 = -0.7872
        du[1] = a*(x[2]-x[1]-(m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))))
        du[2] = x[1]-x[2]+x[3]
        du[3] = -b*x[2]-c*x[3]
    end
    α = [0.93, 0.99, 0.92];
    u0 = [0.2; -0.1; 0.1];
    tspan = (0, 0.5);
    prob = FODEProblem(chua!, α, u0, tspan)
    sol = solve(prob, NonLinearAlg(), dt = 0.1)

    @test isapprox(test_sol(sol), [  0.2   0.11749     0.074388  0.0733938  0.117483   0.210073
    -0.1  -0.0590683  -0.018475  0.0192931  0.0534393  0.0844175
     0.1   0.224134    0.282208  0.286636   0.246248   0.168693]; atol=1e-3)
end
#=
@testset "Test Product Integral Explicit method" begin
    fun(y, p, t)=1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, (0, 5))
    sol=solve(prob, 0.5, PIEX())

    @test isapprox(sol.u, [0.0
    0.5840920370824765
    0.6497603208539517
    0.6965874966593231
    0.7309421360673156
    0.7569810570702951
    0.7773192302138408
    0.7936279277715693
    0.8070043316131614
    0.8181895332990293
    0.8276979913681598]; atol=1e-4)
end
=#



@testset "Test Product Integration explicit method for multi-terms FODE" begin
    tspan = (0, 30); h = 0.01
    rightfun(y, p, x) = 172/125*cos(4/5*x)
    u0 = [0, 0, 0, 0, 0, 0]
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, u0, tspan)
    sol = solve(prob, MTPIEX(), dt=h)

    @test isapprox(sol.u[end-20:end], [0.2441313200142332
    0.24719587414026023
    0.25016860604787894
    0.25304990316172393
    0.2558401698074988
    0.25853982709645607
    0.2611493128113693
    0.2636690812886542
    0.26609960329935234
    0.26844136592362133
    0.2706948724283791
    0.2728606421363793
    0.27493921029876844
    0.2769311279569564
    0.2788369618131128
    0.28065729408669426
    0.28239272237789237
    0.28404385952073596
    0.2856113334413095
    0.28709578700743066
    0.288497877879621]; atol=1e-4)
end

@testset "Test Product integration with predictor-corrector method for multi-terms FODE" begin
    tspan = (0.0, 10.0);  h = 0.5
    rightfun(u, p, x) = 172/125*cos(4/5*x)
    u0 = [0, 0, 0, 0, 0, 0]
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], tspan)
    
    sol = solve(prob, PIPECE(), dt = h)

    @test isapprox(sol.u, [  0.0
    0.01942631486406702
    0.13703330688409168
    0.3665463925810521
    0.6301167829803188
    0.7968019190472511
    0.7407968242020087
    0.3985067557399453
   -0.19598223658593597
   -0.9100047185814064
   -1.5470928325016606
   -1.9063815429019095
   -1.8472382231556892
   -1.3378611193669085
   -0.47174020197161115
    0.5537908253525861
    1.49062009988213
    2.108384867938473
    2.2577394376960034
    1.9096795116394878
    1.160023434854465]; atol=1e-4)
end


@testset "Test implicit Product integration trapezoidal type method for multi-terms FODE with more precise steps" begin
    tspan = (0, 30); h = 0.01
    rightfun(u, p, x) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], tspan)
    sol = solve(prob, PITrap(), dt=h)

    @test isapprox(sol.u[end-20:end], [0.2062614941629048
    0.21034012743472855
    0.21433273142169978
    0.21823954838305873
    0.22206083656949055
    0.22579687013768238
    0.2294479390502145
    0.23301434898224194
    0.23649642122087866
    0.23989449256556789
    0.24320891522328197
    0.24644005670048225
    0.2495882996962482
    0.25265404199064473
    0.2556376963304972
    0.25853969031464596
    0.2613604662731697
    0.2641004811488366
    0.266760206373581
    0.26934012774288163
    0.27184074528802615]; atol=1e-4)
end


@testset "Test Product integration with predictor-corrector method for multi-terms FODE" begin
    tspan = (0, 1); h = 0.01
    rightfun(y, p, x) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], tspan)
    sol = solve(prob, PIPECE(), dt=h)

    @test isapprox(sol.u[end-10:end], [0.12517053205360132
    0.128867333370644
    0.13262214977857692
    0.13643478752690427
    0.14030503437613553
    0.1442326596506723
    0.14821741429463087
    0.15225903093059495
    0.15635722392129503
    0.1605116894342043
    0.16472210550904948]; atol=1e-4)
end



@testset "Test implicit Product integration rectangular type method for multi-terms FODE" begin
    T = 10; h = 0.5
    rightfun(y, p, x) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], (0, T))
    
    sol = solve(prob, PIRect(), dt=h)

    @test isapprox(sol.u, [0.0
    0.01586240520261297
    0.1119084041841651
    0.2954711964704977
    0.5008017075335748
    0.6241431840961683
    0.5697362571789858
    0.2914052658364976
   -0.1836104680656457
   -0.7558230882351726
   -1.278848877489753
   -1.6008246449445735
   -1.6091002088212767
   -1.2648474636427782
   -0.6175431717278957
    0.2046270023806023
    1.0255234402793145
    1.664252313922362
    1.9791244106539438
    1.9009297334987802
    1.447593455242784]; atol=1e-4)
end

@testset "Test implicit Product integration rectangular type method for multi-terms FODE with more precise steps" begin
    T = 30; h = 0.01
    rightfun(x, y) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], (0, T))
    
    sol = solve(prob, PIRect(), dt=h)

    @test isapprox(sol.u[end-20:end], [0.16675726971457058
    0.17176357358985567
    0.17668977521149823
    0.18153598448758948
    0.18630232641110522
    0.19098894098656477
    0.19559598315439394
    0.20012362271686987
    0.2045720442531887
    0.2089414470422814
    0.2132310482116031
    0.21744309867908648
    0.22157681490935055
    0.22563245421407044
    0.22961028818223225
    0.2335106024963484
    0.2373336968672225
    0.24107988495355626
    0.244749494257382
    0.24834286603206457
    0.251860355168466]; atol=1e-4)
end

@testset "Test implicit Product integration trapezoidal type method for multi-terms FODE" begin
    T = 10; h = 0.5
    rightfun(y, p, x) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], (0, T))
    
    sol = solve(prob, PITrap(), dt=h)

    @test isapprox(sol.u, [0.0
    0.02157284502548659
    0.153474720938717
    0.41530591952310447
    0.723333042644268
    0.9282329990906286
    0.879936169829285
    0.4964927100233818
   -0.1888267323884115
   -1.0211143586249927
   -1.7629487467873093
   -2.167707731932248
   -2.062971791641991
   -1.415305444047155
   -0.3530008235484339
    0.86303702756878
    1.9127384461430172
    2.5119503312180154
    2.498678841978473
    1.882454666603612
    0.8406427306579508]; atol=1e-4)
end

@testset "Test product integration PECE method for FODEProblem" begin
    function test(du, u, p, t)
        du[1] = 1-4*u[1]+u[1]^2*u[2]
        du[2] = 3*u[1]-u[1]^2*u[2]
    end
    
    alpha = [0.8, 0.7]
    tspan=(0, 0.1)
    y0=[1.2; 2.8]
    param=[1 3]
    prob = FODEProblem(test, alpha, y0, tspan)
    sol = solve(prob, PECE(), dt=0.01)

    @test isapprox(test_sol(sol), [1.2  1.2061   1.21042  1.21421  1.21767  1.22089  1.22392  1.22678  1.2295   1.2321   1.23459
    2.8  2.78118  2.76953  2.7596   2.75064  2.74235  2.73455  2.72715  2.72008  2.71329  2.70673]; atol=1e-3)
end

@testset "Test explicit product integration method for FODEProblem" begin
    function test(du, u, p, t)
        du[1] = 1-4*u[1]+u[1]^2*u[2]
        du[2] = 3*u[1]-u[1]^2*u[2]
    end
    
    alpha = [0.8, 0.7]
    tspan = (0, 0.2)
    y0=[1.2; 2.8]
    param=[1 3]
    prob = FODEProblem(test, alpha, y0, tspan)
    sol = solve(prob, PIEX(), dt = 0.015)

    @test isapprox(test_sol(sol), [1.2  1.20865  1.21458  1.21975  1.22443  1.22874  1.23276  1.23652  1.24006  1.24339  1.24653  1.2495   1.25229  1.25493  1.25575
    2.8  2.77486  2.75941  2.74621  2.73429  2.72326  2.71291  2.70309  2.69373  2.68477  2.67616  2.66786  2.65986  2.65213  2.64964]; atol=1e-2)
end

@testset "Test NewtonGregory method" begin
    a=1; mu=4
    function brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
    end
    alpha=[0.8; 0.8]
    u0=[0.2; 0.03]
    tspan=(0, 0.5)
    prob = FODEProblem(brusselator, alpha, u0, tspan)
    sol = solve(prob, NewtonGregory(), dt = 0.1)

    @test isapprox(test_sol(sol), [ 0.2   0.200531  0.201161  0.201809  0.202453  0.203089
    0.03  0.165554  0.265678  0.355635  0.439544  0.519211]; atol=1e-4)
end

@testset "Test FLMMBDF" begin
    a=1; mu=4
    function brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
    end
    alpha=[0.8; 0.8]
    u0=[0.2; 0.03]
    tspan=(0, 0.5)
    prob = FODEProblem(brusselator, alpha, u0, tspan)
    sol = solve(prob, BDF(), dt = 0.1)

    @test isapprox(test_sol(sol), [ 0.2   0.200531  0.201161  0.201809  0.202453  0.203089
    0.03  0.165554  0.265678  0.355635  0.439544  0.519211]; atol=1e-4)
end

@testset "Test Trapezoid method" begin
    a=1; mu=4
    function brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
    end
    alpha=[0.8; 0.8]
    u0=[0.2; 0.03]
    tspan=(0, 0.5)
    prob = FODEProblem(brusselator, alpha, u0, tspan)
    sol = solve(prob, Trapezoid(), dt = 0.1)

    @test isapprox(test_sol(sol), [0.2   0.200531  0.201161  0.201808  0.202452  0.203088
    0.03  0.165554  0.265678  0.355635  0.439545  0.519211]; atol=1e-4)
end


@testset "Test Newton Polynomial" begin
    t0=0; tfinal=5
    tspan = (0.0, 5.0)
    α = [0.98, 0.98, 0.98]
    u0 = [-1.0, 1.0, 1.0]
    function fun(du, u, p, t)
        b=0.1;
        du[1] = -sin(u[2])-b*u[1]
        du[2] = -sin(u[3])-b*u[2]
        du[3] = -sin(u[1])-b*u[3]
    end
    prob = FODEProblem(fun, α, u0, tspan)
    sol = solve(prob, NewtonPolynomial(), dt = 0.5)

    @test isapprox(test_sol(sol), [-1.0  -1.37074   -1.46124    -1.21771   -0.884241  -0.49327   -0.0897447   0.338635   0.793532   1.2323    1.55187
    1.0   0.529265  -0.0101035  -0.430956  -0.733314  -0.916321  -1.06158    -1.2647    -1.5767    -1.99613  -2.37772
    1.0   1.37074    1.8176      2.17624    2.48899    2.67047    2.6624      2.46322    2.07376    1.55057   1.0049]; atol=1e-4)
end

#=
@testset "Test Atangana Seda method" begin
    fun(y, p, t) = 1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, (0, 5))
    sol = solve(prob, 0.5, AtanganaSeda())

    @test isapprox(sol.u, [0.0
    0.5
    0.625
    0.4122403564148137
    0.9011748394667164
   -0.11601114237053156
    2.3113598236715522
   -3.0656444401780365
    9.199001384099954
  -18.488162171353544
   44.23628311938973]; atol=1e-3)
end
=#
@testset "Test Atangana Seda method for Atangana-Baleanu FODEProblem" begin
    tspan = (0.0, 2.0)
    h=0.5
    α = [0.99, 0.99, 0.99]
    u0 = [-0.1; 0.1; -0.1]
    function fun(du, u, p, t)
        a=0.2;b=4;ζ=8;δ=1;
        du[1] = -u[2]^2-u[3]^2-a*u[1]+a*ζ
        du[2] = u[1]*u[2]-b*u[1]*u[3]-u[2]+δ
        du[3] = b*u[1]*u[2]+u[1]*u[3]-u[3]
    end
    prob = FODEProblem(fun, α, u0, tspan)
    sol = solve(prob, AtanganaSedaAB(), dt = 0.5)

    @test isapprox(test_sol(sol), [ -0.1   0.7    1.18511  -1.41522  -31.0872
    0.1   0.525  1.08088  -4.03383   32.0085
   -0.1  -0.065  1.03463   4.10537   13.3441]; atol=1e-2)
end

@testset "Test AtanganaSedaCF method" begin
    a₁, a₂, a₃, a₄, a₅, a₆, a₇ = 3, 0.5, 4, 3, 4, 9, 4
    function LotkaVolterra(du, u, p, t)
        du[1] = u[1]*(a₁-a₂*u[1]-u[2]-u[3])
        du[2] = u[2]*(1-a₃+a₄*u[1])
        du[3] = u[3]*(1-a₅+a₆*u[1]+a₇*u[2])
    end
    u0 = [0.5, 0.9, 0.1]
    tspan = (0.0, 2.0)
    α = ones(3)*0.98;
    prob = FODEProblem(LotkaVolterra, α, u0, tspan)
    sol = solve(prob, AtanganaSedaCF(), dt = 0.5)

    isapprox(sol.u, [ 0.0   1.375   3.37602   -5.46045    449.712
    0.0  -0.45   -0.492188  -3.4828      62.7789
    0.0   0.61    3.94806   96.048    -5989.62]; atol=1e-3)
end