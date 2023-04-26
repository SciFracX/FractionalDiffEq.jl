# Commensurate case

using FractionalDiffEq, Plots

function RF(du, u, t)
    du[1] = u[2]*(u[3]-1+u[1]*u[1])+0.1*u[1]
    du[2] = u[1]*(3*u[3]+1-u[1]*u[1])+0.1*u[2]
    du[3] = -2*u[3]*(0.98+u[1]*u[2])
end

LE=FOLyapunov(RF, [0.999, 0.999, 0.999], 0, 0.02, 300, [0.1; 0.1; 0.1], 0.005, 1000)
plot(LE)
plot(tspan, LE[1, :])
plot!(tspan, LE[2, :])
plot!(tspan, LE[3, :])



# Noncommensurate case

using FractionalDiffEq, Plots

function LE_RF_TEST(du, u, p, t)
    du[1] = u[2]*(u[3]-1+u[1]^2) + 0.1*u[1]
    du[2] = u[1]*(3*u[3]+1-u[1]^2) + 0.1*u[2]
    du[3] = -2*u[3]*(0.98+u[1]*u[2])
end
LE = FOLyapunov(LE_RF_TEST, [0.995, 0.992, 0.996], 0, 0.1, 1000, [1,1,1], 0.01, 1000)
plot(LE)