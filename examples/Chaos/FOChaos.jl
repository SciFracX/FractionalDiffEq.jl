# Fractional Order Chaotic examples come from 
#=
@inproceedings{Petr2011FractionalOrderNS,
title={Fractional-Order Nonlinear Systems: Modeling, Analysis and Simulation},
author={Ivo Petr{\'a}{\vs}},
year={2011}
}
=#
#

#=======Fractional Order Liu System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [0.2, 0, 0.5]
tf=100
function Liu!(du, u, p, t)
    a, b, c, e, n, m = 1, 2.5, 5, 1, 4, 4
    du[1] = -a*u[1]-e*u[2]^2
    du[2] = b*u[2]-n*u[1]*u[3]
    du[3] = -c*u[3]+m*u[1]*u[2]
end
prob = FODESystem(Liu!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Liu System")


#=======Fractional Order Duffing System=======#
using FractionalDiffEq

h=0.005
alpha = [0.9, 1]
x0 = [0.21, 0.31]
tf=100
function Duffing!(du, u, p, t)
    α, δ, ω = 0.15, 0.3, 1
    du[1] = u[2]
    du[2] = u[1]-α*u[2]-u[1]^3+δ*cos(ω*t)
end
prob = FODESystem(Duffing!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Duffing System")


#=======Fractional Order Chen System=======#
using FractionalDiffEq

h=0.005
alpha = [0.9, 0.9, 0.9]
x0 = [-9, -5, 14]
tf=100
function Chen!(du, u, p, t)
    a, b, c, d = 35, 3, 28, -7
    du[1] = a*(u[2]-u[1])
    du[2] = d*u[1]-u[1]*u[3]+c*u[2]
    du[3] = u[1]*u[2]-b*u[3]
end
prob = FODESystem(Chen!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Chen System")


#=======Fractional Order Van der Pol Oscillator=======#
using FractionalDiffEq

h=0.005
alpha = [1.2, 0.8]
x0 = [0.2, -0.2]
tf=60
function VanderPol!(du, u, p, t)
    ϵ = 1
    du[1] = u[2]
    du[2] = -u[1]-ϵ*(u[1]^2-1)*u[2]
end
prob = FODESystem(VanderPol!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Van der Pol Oscillator")


#=======Fractional Order Lotka Volterra System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [1, 1.4, 1]
tf=50
function LotkaVolterra!(du, u, p, t)
    a, b, c, d, e, p, s = 1, 1, 1, 1, 2, 3, 2.7
    du[1] = a*u[1]+e*u[1]^2-b*u[1]*u[2]-s*u[3]*u[1]^2
    du[2] = -c*u[2]+d*u[1]*u[2]
    du[3] = -p*u[3]+s*u[3]*u[1]^2
end
prob = FODESystem(LotkaVolterra!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lotka Volterra System")


#=======Fractional Order Lu System=======#
using FractionalDiffEq

h=0.005
alpha = [0.985, 0.99, 0.98]
x0 = [0.2, 0.5, 0.3]
tf=60
function Lu!(du, u, p, t)
    a, b, c = 36, 3, 20
    du[1] = a*(u[2]-u[1])
    du[2] = -u[1]*u[3]+c*u[2]
    du[3] = u[1]*u[2]-b*u[3]
end
prob = FODESystem(Lu!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lu System")


#=======Fractional Order Rossler System=======#
using FractionalDiffEq

h=0.005
alpha = [0.9, 0.85, 0.95]
x0 = [0.5, 1.5, 0.1]
tf=120
function Rossler!(du, u, p, t)
    a, b, c = 0.5, 0.2, 10
    du[1] = -u[2]-u[3]
    du[2] = u[1]+a*u[2]
    du[3] = b+u[3]*(u[1]-c)
end
prob = FODESystem(Rossler!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Rossler System")


#=======Fractional Order Arneodo System=======#
using FractionalDiffEq

h=0.005
alpha = [0.97, 0.97, 0.96]
x0 = [-0.2, 0.5, 0.2]
tf=100
function Arneodo!(du, u, p, t)
    b1, b2, b3, b4 = -5.5, 3.5, 0.8, -1.0
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -b1*u[1]-b2*u[2]-b3*u[3]+b4*u[1]^3
end
prob = FODESystem(Arneodo!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Arneodo System")


#=======Fractional Order Genesio-Tesi System=======#
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [-0.1, 0.5, 0.2]
tf=100
function GenesioTesi!(du, u, p, t)
    b1, b2, b3, b4 = 1, 1.1, 0.4, 1.0
    du[1] = u[2]
    du[2] = u[3]
    du[3] -b1*u[1]-b2*u[2]-b3*u[3]+b4*u[1]^2
end
prob = FODESystem(GenesioTesi!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Genesio-Tesi System")