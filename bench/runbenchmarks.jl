using BenchmarkTools: BenchmarkTools, BenchmarkGroup, @btime, @benchmarkable
using FractionalDiffEq
using Statistics: median

@info sprint(versioninfo)

const SUITE = BenchmarkGroup()

α = [0.8, 0.8]
u0 = [0.2, 0.03]

function Brusselator!(du, u, p, t)
    a, μ = 1, 4
    du[1] = a-(μ+1)*u[1]+(u[1])^2*u[2]
    du[2] = μ*u[1]-(u[1])^2*u[1]
end

prob = FODESystem(Brusselator!, α, u0, (0, 20))

SUITE["FLMM"]["BDF"] = @benchmarkable solve(prob, BDF(), dt=0.01)
SUITE["FLMM"]["Trapezoid"] = @benchmarkable solve(prob, Trapezoid(), dt=0.01)
SUITE["FLMM"]["NewtonGregory"] = @benchmarkable solve(prob, NewtonGregory(), dt=0.01)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

BenchmarkTools.save(joinpath(@__DIR__, "benchmark_results.json"), median(results))