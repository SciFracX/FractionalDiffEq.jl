@testset "Test Diethelm PECE algorithms" begin
    fun(x, y) = 1-y
    prob = SingleTermFODEProblem(fun, 1.8, 0, 5)
    result = solve(prob, 0.01, PECE())
    tspan = collect(0:0.01:5)

    target = []

    for i in 0:0.01:5
        push!(target, i^1.8*mittleff(1.8, 2.8,-i^1.8))
    end

    @test isapprox(result.u, target; atol=1)
end

@testset "Test Matrix discrete method" begin
    fun(x, y) = 1-y
    u0 = [0, 0]
    tspan=collect(0.01:0.01:5)
    target = []
    for i in 0.01:0.01:5
        push!(target, i^0.5*mittleff(0.5, 1.5,-i^0.5))
    end

    singleprob = MultiTermsFODEProblem([1, 1], [0.5, 0], 1, u0, 5)
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

    highsingleprob = MultiTermsFODEProblem([1, 1], [1.8, 0], 1, u0, 20)
    yaresult = solve(highsingleprob, 0.01, FODEMatrixDiscrete())

    @test isapprox(yaresult.u, yatarget; atol=1)

end

@testset "Test Closed Form method" begin
    t=collect(0:0.002:10);
    u0 = [0, 0, 0, 0]
    rightfun(x) = sin(x)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3], u0, t)
    result=solve(prob, ClosedForm())
end

@testset "Test ClosedFormHankelM method" begin
    t = collect(0:0.5:1);
    u0 = [0, 0, 0, 0]
    rightfun(x)=sin(x^2)
    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3], u0, t)
    result = solve(prob, ClosedFormHankelM())

    @test result.u ≈ [0.0; 0.08402140107687359; 0.3754974742112727]
end

@testset "Test GL method for FODESystem" begin
    h=0.5; tf=1
    alpha = [0.99, 0.99, 0.99]
    x0 = [1, 0, 1]
    function testf!(du, u, p, t)
        a, b, c = 10, 28, 8/3
        du[1] = a*(u[2]-u[1])
        du[2] = u[1]*(b-u[3])-u[2]
        du[3] = u[1]*u[2]-c*u[3]
    end
    prob = FODESystem(testf!, alpha, x0, tf)
    result = solve(prob, h, GL())
    @test isapprox(result, [1.0 0.0 1.0
    -4.04478 13.5939 -0.352607
    84.8074 -51.1251 -27.5541]; atol=1e-4)
end

@testset "Test GL method" begin
    fun(x, y) = 1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 1)
    sol = solve(prob, 0.1, GL())

    @test isapprox(sol.u, [0.0
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
    x0 = [0.2; -0.1; 0.1];
    h = 0.1; tn = 0.5;
    prob = FODESystem(chua!, α, x0, tn)
    result = solve(prob, h, NonLinearAlg())

    @test isapprox(result, [ 0.2 -0.1 0.1
    0.11749    -0.0590683  0.224134
    0.074388   -0.018475   0.282208
    0.0733938   0.0192931  0.286636
    0.117483    0.0534393  0.246248
    0.210073    0.0844175  0.168693]; atol=1e-3)
end

@testset "Test Product Integral Explicit method" begin
    fun(t, y)=1-y
    prob = SingleTermFODEProblem(fun, 0.5, 0, 5)
    sol=solve(prob, 0.5, PIEX())

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
    sol=solve(prob, 0.1, PIIM())

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

@testset "Test ChebSpectral method" begin
    n = 10
    rightfun(x) = gamma(7+9/17)/gamma(7+9/17-0.5)*(1+x)^(6+9/17-0.5)
    prob = SingleTermFODEProblem(rightfun, 0.5, 0, 1)
    sol = solve(prob, n, ChebSpectral())

    @test isapprox(sol.u, [-2.103832955612273e-8
    7.611080159690798e-5
    0.01082543428432764
    0.2878301921177669
    2.8447345612767694
   14.117972832745936
   40.9991958543131
   75.6356722636105
   92.37379693143299]; atol=1e-4)
end

@testset "Test Product Integration explicit method for multi-terms FODE" begin
    #T = 10; h = 0.5
    rightfunex(x, y) = 172/125*cos(4/5*x)
    PIEXMultiprob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfunex, [0, 0, 0, 0, 0, 0], 0, 10)
    PIEXMultisol = solve(PIEXMultiprob, 0.5, PIEX())

    @test isapprox(PIEXMultisol.u, [  0.0
    0.028666666666666663
    0.20802154481466284
    0.5777512805039694
    1.032617182281499
    1.3542099549416406
    1.2992204949109856
    0.7162691469691946
   -0.35553236137048394
   -1.6531692565209024
   -2.748115186408473
   -3.188228316718636
   -2.6777963530465527
   -1.226608005191384
    0.7994418493409547
    2.76345190855897
    3.9650094859689915
    3.9040343901844707
    2.5017274305718367
    0.18249065149394283
   -2.23584183347657]; atol=1e-4)
end

@testset "Test Product integration with predictor-corrector method for multi-terms FODE" begin
    T = 10; h = 0.5
    rightfun(x, y) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], 0, T)
    
    sol = solve(prob, h, PIPECE())

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


@testset "Test implicit Product integration rectangular type method for multi-terms FODE" begin
    T = 10; h = 0.5
    rightfun(x, y) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], 0, T)
    
    sol = solve(prob, h, PIIMRect())

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

@testset "Test implicit Product integration trapezoidal type method for multi-terms FODE" begin
    T = 10; h = 0.5
    rightfun(x, y) = 172/125*cos(4/5*x)
    prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], 0, T)
    
    sol = solve(prob, h, PIIMTrap())

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

@testset "Test product integration PECE method for FODESystem" begin
    function testf(t, y)
        p=[1 3]
        return [p[1]-(p[2]+1)*y[1] + y[1]^2*y[2]; p[2]*y[1]-y[1]^2*y[2]]
    end
    
    alpha = [0.8, 0.7]
    t0=0; T=0.1; h=0.01
    y0=[1.2; 2.8]
    param=[1 3]
    prob = FODESystem(testf, alpha, y0, t0, T)
    (t, y) = solve(prob, h, PECE())

    @test isapprox(y, [1.2  1.2061   1.21042  1.21421  1.21767  1.22089  1.22392  1.22678  1.2295   1.2321   1.23459
    2.8  2.78118  2.76953  2.7596   2.75064  2.74235  2.73455  2.72715  2.72008  2.71329  2.70673]; atol=1e-3)
end

@testset "Test explicit product integration method for FODESystem" begin
    function testf(t, y)
        p=[1 3]
        return [p[1]-(p[2]+1)*y[1] + y[1]^2*y[2]; p[2]*y[1]-y[1]^2*y[2]]
    end
    
    alpha = [0.8, 0.7]
    t0=0; T=0.1; h=0.01
    y0=[1.2; 2.8]
    param=[1 3]
    prob = FODESystem(testf, alpha, y0, t0, T)
    (t, y) = solve(prob, h, PIEX())

    @test isapprox(y, [ 1.2  1.20626  1.21061  1.21444  1.21793  1.22118  1.22424  1.22713  1.22989  1.23252  1.23504
    2.8  2.78107  2.76943  2.75949  2.75052  2.74221  2.7344   2.72697  2.71988  2.71305  2.70647]; atol=1e-4)
end

@testset "Test FLMMNewtonGregory method" begin
    a=1; mu=4
    function Brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
        du
    end
    alpha=[0.8; 0.8]
    t0=0; tfinal=0.5; y0=[0.2; 0.03]
    h=0.1
    prob = FODESystem(Brusselator, alpha, y0, t0, tfinal)
    (t, y) = solve(prob, h, FLMMNewtonGregory())

    @test isapprox(y, [ 0.2   0.200531  0.201161  0.201809  0.202453  0.203089
    0.03  0.165554  0.265678  0.355635  0.439544  0.519211]; atol=1e-4)
end

@testset "Test FLMMBDF" begin
    a=1; mu=4
    function Brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
        du
    end
    alpha=[0.8; 0.8]
    t0=0; tfinal=0.5; y0=[0.2; 0.03]
    h=0.1
    prob = FODESystem(Brusselator, alpha, y0, t0, tfinal)
    (t, y) = solve(prob, h, FLMMBDF())

    @test isapprox(y, [ 0.2   0.200531  0.201161  0.201809  0.202453  0.203089
    0.03  0.165554  0.265678  0.355635  0.439544  0.519211]; atol=1e-4)
end

@testset "Test FLMMTrap" begin
    a=1; mu=4
    function Brusselator(du, u, p, t)
        du[1] = a-(mu+1)*u[1]+u[1]^2*u[2]
        du[2] = mu*u[1]-u[1]^2*u[2]
        du
    end
    alpha=[0.8; 0.8]
    t0=0; tfinal=0.5; y0=[0.2; 0.03]
    h=0.1
    prob = FODESystem(Brusselator, alpha, y0, t0, tfinal)
    (t, y) = solve(prob, h, FLMMTrap())

    @test isapprox(y, [0.2   0.200531  0.201161  0.201808  0.202452  0.203088
    0.03  0.165554  0.265678  0.355635  0.439545  0.519211]; atol=1e-4)
end