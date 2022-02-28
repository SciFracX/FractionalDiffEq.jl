using FractionalDiffEq, Plots, SpecialFunctions

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

result = solve(3, 0.5, 1, 20, 1, 20, 1, 1, fx0, fgz, f0t, flt, ADV_DIF())
XX, YY = meshgrid(0.05^2/6 .*(0:21), 0:0.05:1)
plotlyjs()
plot(XX, YY, result, st=:surface)