using FractionalDiffEq, Plots
T = 30; h = 0.01
rightfun(x, y) = 172/125*cos(4/5*x)
prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], 0, T)
sol = solve(prob, h, PIEX())
plot(sol, legend=:bottomright)


    T = 10; h = 0.5
    rightfunex(x, y) = 172/125*cos(4/5*x)
    PIEXMultiprob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfunex, [0, 0, 0, 0, 0, 0], 0, 10)
    PIEXMultisol = solve(PIEXMultiprob, 0.5, PIEX())