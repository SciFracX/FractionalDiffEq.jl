using FractionalDiffEq
using Test

@testset "Test Lyapunov exponents" begin
    function RF(du, u, t)
        du[1] = u[2]*(u[3]-1+u[1]*u[1])+0.1*u[1]
        du[2] = u[1]*(3*u[3]+1-u[1]*u[1])+0.1*u[2]
        du[3] = -2*u[3]*(0.98+u[1]*u[2])
    end
    
    (LE, tspan)=FOLyapunov(RF, 0.999, 0, 0.02, 300, [0.1; 0.1; 0.1], 0.005, 1000)

    # Test the latest computed(final) Lyapunov exponents
    @test isapprox(LE[end-2:end], [0.06111650166568285
    0.0038981396237095034
   -1.8324646820425692]; atol=1e-3)
end
