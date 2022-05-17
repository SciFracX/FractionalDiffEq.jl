# An example of system of Fractional Difference Equations

```math
\begin{cases}
{_3^G\nabla}_k^\alpha x_1(k+1)=-0.05x_2(k)-0.05x_3(k)+0.01\tanh(x_2(k))\\
{_3^G\nabla}_k^\alpha x_2(k+1)=0.05x_1(k)+0.02x_2(k)+0.01\tanh(x_1(k))\\
{_3^G\nabla}_k^\alpha x_3(k+1)=0.1-0.2x_3(k)+0.05x_1(k)x_3(k)+0.01\tanh(x_3(k))
\end{cases}
```


```julia
using FractionalDiffEq, Plots

function sys!(du, u, p, t)
    du[1] = -0.05*u[2] - 0.05*u[3] + 0.01*tanh(u[2])
    du[2] = 0.05*u[1] + 0.02*u[2] + 0.01*tanh(u[1])
    du[3] = 0.1 - 0.2*u[3] + 0.05*u[1]*u[3] + 0.01*tanh(u[3])
end
prob = FractionalDifferenceSystem(sys!, 0.98, [1, -1, 0])
result = solve(prob, 7, GL())

plot(result[1, :], result[2, :], result[3, :], seriestype=:scatter)
```

![FNetwork](./assets/fractionalneuralnetwork.png)