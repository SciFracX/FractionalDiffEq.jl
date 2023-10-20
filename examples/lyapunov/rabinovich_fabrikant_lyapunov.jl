# Commensurate case

using FractionalDiffEq, Plots

function RF(du, u, p, t)
    a, b, c = p
    du[1] = u[2]*(u[3]-a+u[1]*u[1])+b*u[1]
    du[2] = u[1]*(3*u[3]+1-u[1]*u[1])+b*u[2]
    du[3] = -2*u[3]*(c+u[1]*u[2])
end

p = [1, 0.1, 0.98]

LE=FOLyapunov(RF, [0.999, 0.999, 0.999], 0, 0.02, 300, [0.1; 0.1; 0.1], 0.005, 1000, p)
plot(LE)

# Noncommensurate case

using FractionalDiffEq, Plots

function LE_RF_TEST(du, u, p, t)
    du[1] = u[2]*(u[3]-1+u[1]^2) + 0.1*u[1]
    du[2] = u[1]*(3*u[3]+1-u[1]^2) + 0.1*u[2]
    du[3] = -2*u[3]*(0.98+u[1]*u[2])
end
LE = FOLyapunov(LE_RF_TEST, [0.995, 0.992, 0.996], 0, 0.1, 1000, [1,1,1], 0.01, 1000)
plot(LE)

