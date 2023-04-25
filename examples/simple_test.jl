using ForwardDiff, FractionalDiffEq

function FO_NC_Lyapunov(fun, order, t_start, h_norm, t_end, u0, h, out)# TODO: Generate the Lyapunov exponent plot
    ne::Int = length(u0) # System dimension

    tspan = Float64[]
    LE = Float64[]

    # Generate extend system with jacobian
    function extend_fun(du, u, p, t)
        current_du = zeros(ne)
        fun(current_du, u[1:ne], p, t)
        extend_du = jacobian_of_fdefun(fun, t, u[1:ne], p)
        extend_du_mul_u = extend_du*reshape(u[ne+1:end], ne, ne)
        du[:] = [current_du[:]; extend_du_mul_u[:]]
    end

    x = zeros(Float64, ne*(ne+1))
    x0 = zeros(Float64, ne*(ne+1))
    c = zeros(Float64, ne)
    gsc = zeros(Float64, ne)
    zn = zeros(Float64, ne)
    n_it = round(Int, (t_end-t_start)/h_norm)
    x[1:ne] = u0
    q = repeat(order, ne+1, 1) # fractional order of the extend system
    for i=1:ne
        x[(ne+1)*i]=1.0
    end
    t = t_start
    LExp = zeros(ne)
    for it=1:n_it
        prob = FODESystem(extend_fun, q[:], x[:], (t, t+h_norm))
        sol = solve(prob, h, PECE())
        Y = sol.u
        t = t+h_norm
        Y = Y'

        x = Y[size(Y, 1), :]
        for i=1:ne
            for j=1:ne
                x0[ne*i+j]=x[ne*j+i]
            end
        end
        zn[1] = 0.0
        for j=1:ne
            zn[1] = zn[1]+x0[ne*j+1]^2
        end
        zn[1] = sqrt(zn[1])
        for j=1:ne
            x0[ne*j+1] = x0[ne*j+1]/zn[1]
        end
        for j=2:ne
            for k=1:(j-1)
                gsc[k] = 0.0
                for l=1:ne
                    gsc[k] = gsc[k]+x0[ne*l+j]*x0[ne*l+k]
                end
            end
            for k=1:ne
                for l=1:j-1
                    x0[ne*k+j]=x0[ne*k+j]-gsc[l]*x0[ne*k+l]
                end
            end
            zn[j]=0.0
            for k=1:ne
                zn[j]=zn[j]+x0[ne*k+j]^2
            end
            zn[j]=sqrt(zn[j])
            for k=1:ne
                x0[ne*k+j]=x0[ne*k+j]/zn[j]
            end
        end
        for k=1:ne
            c[k]=c[k]+log(zn[k])
            LExp[k]=c[k]/(t-t_start)
        end
        for i=1:ne
            for j=1:ne
                x[ne*j+i]=x0[ne*i+j]
            end
        end
        x=transpose(x)
        mod(it, out)==0 ? println(LExp) : nothing
        LE = [LE; LExp]
        tspan = [tspan; t]
    end
    LE = reshape(LE, ne, :)
    return LE, tspan
end

function jacobian_of_fdefun(f, t, y, p)
    ForwardDiff.jacobian(y) do y
    du = similar(y)
    f(du, y, p, t)
    du
    end
end

function LE_RF_TEST(du, u, p, t)
    du[1] = u[2]*(u[3]-1+u[1]^2) + 0.1*u[1]
    du[2] = u[1]*(3*u[3]+1-u[1]^2) + 0.1*u[2]
    du[3] = -2*u[3]*(0.98+u[1]*u[2])
end
LE, tspan = FO_NC_Lyapunov(LE_RF_TEST, [0.995, 0.992, 0.996], 0, 0.1, 1000, [1,1,1], 0.01, 1000)
using Plots

plot(tspan, LE[1, :])
plot!(tspan, LE[2, :])
plot!(tspan, LE[3, :])
savefig("D:\\Project\\noncommensurate_lyapunov.png")