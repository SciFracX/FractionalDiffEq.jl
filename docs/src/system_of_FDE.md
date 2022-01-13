# System of fractional differential equations

Many "real life" situations are governed by a system of fractional differential equations.

So here, we will look at an example: Chua circuit.

The circuit diagram of the Chua system is shown below:

![Chua diode](./assets/chua_diode.svg)

> Here, **``N_R``** is the [memoristor](https://en.wikipedia.org/wiki/Memristor), which is a non-linear electrical component relating electric charge and magnetic flux linkage.

Let's see if we abstract the Chua system into a fractional differential equation system:

```math
\begin{cases}
D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x+1|-|x-1|)]\\
D^{\alpha_2}y=x-y+z\\
D^{\alpha_3}z=-10.593y-0.268z
\end{cases}
```

Use the ```NonLinearAlg``` algorithm in FractionalDiffEq.jl to solve the Chua system and plot the result:

```julia
using FractionalDiffEq
using Plots

function chua(t, x, k)
    a = 10.725
    b = 10.593
    c = 0.268
    m0 = -1.1726
    m1 = -0.7872

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
h = 0.01;
tn = 200;
result = solve(chua, alpha, x0, h, tn, NonLinearAlg())

gr()
plot(result[:, 1], result[:, 2], title="Chua System", legend=:bottomright)
```

![Chua](./assets/chua.png)

Cheers!🎉🎉🎉

It is noteworthy that in the reference book Fractional Calculus and Fractional-order Control[^1], the computing time is almost 20 minutes to solve this problem in [FOTF toolbox](https://www.mathworks.com/matlabcentral/fileexchange/60874-fotf-toolbox), while in FractionalDiffEq.jl, the computing time has a speed up of about 2 times, only cost 8 minutes and 31 seconds!!

[^1]: 分数阶微积分学与分数阶控制 薛定宇 ISBN:9787030543981 Page 208