using FractionalDiffEq
using Plots

function chua(t, x, k)
    a=10.725
    b=10.593
    c=0.268
    m0=-1.1726
    m1=-0.7872

    if k == 1
        f = m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))
        y = a*(x[2]-x[1]-f)
        return y
    elseif k == 2
        y = x[1]-x[2]+x[3]
        return y
    elseif k == 3
        y = -b*x[2]-c*x[3]
        return y
    end
end

alpha = [0.93, 0.99, 0.92];
x0 = [0.2; -0.1; 0.1];
h = 0.001;
tn = 500;
result = solve(chua, alpha, x0, h, tn, NonLinearAlg(), 10000)

gr()
plot(result[:, 1], result[:, 2], title="Chua System with Short Memory Effect", legend=:bottomright)