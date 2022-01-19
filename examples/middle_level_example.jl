using FractionalDiffEq, Plots, SpecialFunctions

α = 0.5

function testfun(t, y)
    return 40320/gamma(9-α)*t^(8-α) - 3*gamma(5+α/2)/gamma(5-α/2)*t^(4-α/2) + 9/4*gamma(α+1) + (3/2*t^α/2-t^4)^3-y^1.5
end

prob = FODEProblem(testfun, α, 0, 1, 0.001)

result = solve(prob, PECE())

tfun(t) = t^8-3t^4.25+9/4*t^0.5

tspan = collect(0:0.001:1)

plot(tspan, result)

plot!(tspan, tfun, lw=2)