using FractionalDiffEq, Plots
plotlyjs()
(X, Y) = meshgrid(collect(-6:0.2:6), collect(-6:0.2:6))

Z = X.+Y.*im

L = mittleff(0.8, 0.9, Z)
L1=real.(L)
id = findall(x -> x>3, L1)
L1[id].=3
ix = findall(x-> x<-3, L1)
L1[ix].=-3
plot(X, Y, L1, st=:surface)