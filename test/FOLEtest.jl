using FractionalDiffEq
using Test

@testset "Test Commensurate Lyapunov exponents" begin
    function RF(du, u, t)
        du[1] = u[2]*(u[3]-1+u[1]*u[1])+0.1*u[1]
        du[2] = u[1]*(3*u[3]+1-u[1]*u[1])+0.1*u[2]
        du[3] = -2*u[3]*(0.98+u[1]*u[2])
    end

    LE = FOLyapunov(RF, [0.999, 0.999, 0.999], 0, 0.02, 300, [0.1; 0.1; 0.1], 0.005, 1000)

    # Test the latest computed(final) Lyapunov exponents
    @test isapprox(LE.LE[end-2:end], [0.06111650166568285
    0.0038981396237095034
   -1.8324646820425692]; atol=1e-3)
end

@testset "Test Noncommensurate Lyapunov exponents" begin
    function LE_RF_TEST(du, u, p, t)
        du[1] = u[2]*(u[3]-1+u[1]^2) + 0.1*u[1]
        du[2] = u[1]*(3*u[3]+1-u[1]^2) + 0.1*u[2]
        du[3] = -2*u[3]*(0.98+u[1]*u[2])
    end

    LE = FOLyapunov(LE_RF_TEST, [0.995, 0.992, 0.996], 0, 0.1, 1000, [1,1,1], 0.01, 1000)   
    @test isapprox(LE.LE[end-2:end], [-0.0033493403496077058
    -0.0050780674255360226
    -1.774132087623718]; atol=1e-3)
end