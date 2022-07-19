using FractionalDiffEq, Plots

τ=3.1416; h=0.01; α=0.4; tspan = (0, 70)
y0(t) = [sin(t)*cos(t); sin(t)*cos(t); cos(t)^2-sin(t)^2; cos(t)^2-sin(t)^2]
A=[0 0 1 0; 0 0 0 1; 0 -2 0 0; -2 0 0 0]
B=[0 0 0 0; 0 0 0 0; -2 0 0 0; 0 -2 0 0]
f=[0; 0; 0; 0]

prob = FDDEMatrixProblem(α, τ, A, B, f, y0, tspan)
sol=solve(prob, h, MatrixForm())
plot(sol[:, 1], sol[:, 3])