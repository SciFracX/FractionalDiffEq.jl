using FractionalDiffEq, Plots, SpecialFunctions

phi(x) = 0
α = 0.9;
τ = 0.1;
function f(t, y, ϕ)
    2 / gamma(3 - α) * t^(2 - α) - 1 / gamma(2 - α) * t^(1 - α) + 2 * τ * t - τ^2 - τ - y +
    ϕ
end
h = 1e-2
realfun(t) = t^2 - t

prob = FDDEProblem(f, phi, α, τ, (0, 5))
sol1 = solve(prob, h, DelayPI())
sol2 = solve(prob, h, DelayABM())
sol3 = sove(prob, h, DelayPECE())
