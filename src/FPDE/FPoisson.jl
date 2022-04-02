mutable struct Mesh
    t
    p
    e
end


function rect_grid2(xmin, xmax, ymin, ymax, nx, ny)
    nt = 2*nx*ny
    np = (nx+1)*(ny+1)


    mesh = Mesh(0, 0, 0)
    mesh.t = zeros(nt, 3)
    mesh.p = zeros(np, 2)
    
    
    nxp1 = nx + 1
    nyp1 = ny + 1
    

    nt  = 0
    for ix = 1:nx
       for iy = 1:ny
          
          iv  = (ix-1)*nyp1 + iy
          iv1 = iv + nyp1
          
          nt = nt + 1
          mesh.t[nt, 1] = iv
          mesh.t[nt, 2] = iv1
          mesh.t[nt, 3] = iv1+1
    
          nt = nt + 1
          mesh.t[nt, 1] = iv
          mesh.t[nt, 2] = iv1+1
          mesh.t[nt, 3] = iv+1
      end
    end
    
    
    
    
    hx   = (xmax-xmin)/nx
    hy   = (ymax-ymin)/ny
    x    = xmin
    
    for ix = 1:nx
    

      i1 = (ix-1)*(ny+1)+1
      i2 = ix*(ny+1)
      mesh.p[i1:i2, 1] = x*ones(nyp1, 1)
      mesh.p[i1:i2, 2] = collect(ymin:hy:ymax)
       
      x = x + hx
    end
    

    i1 = nx*(ny+1)+1
    i2 = (nx+1)*(ny+1)
    mesh.p[i1:i2, 1] = xmax*ones(nyp1,1)
    mesh.p[i1:i2, 2] = collect(ymin:hy:ymax)

    
    mesh.e = ones(2*(nx+ny), 3)
    
    mesh.e[1:ny, 1] = collect(1:ny)
    mesh.e[1:ny, 2] = collect(2:ny+1)
    
    mesh.e[ny+1:nx+ny, 1] = collect(ny+1:ny+1:np-1)
    mesh.e[ny+1:nx+ny, 2] = collect(2*(ny+1):ny+1:np)
    

    mesh.e[nx+ny+1:nx+2*ny, 1] = collect(np-ny:np-1)
    mesh.e[nx+ny+1:nx+2*ny, 2] = collect(np-ny+1:np)

    mesh.e[nx+2*ny+1:2*(nx+ny), 1] = collect(1:ny+1:np-2*ny-1)
    mesh.e[nx+2*ny+1:2*(nx+ny), 2] = collect(ny+2:ny+1:np-ny)
    return mesh
end


function Stiff_Mass(mesh)

    mesh.e[:, 3] .= 1

    N::Int = maximum(maximum(mesh.t[:, 1:3]))


    dirichlet = mesh.e[findall(x->x==1, mesh.e[:, 3]), 1:2]
    dirichlet = unique(sort(dirichlet[:]))
    FreeNodes = setdiff(1:N, dirichlet)



    (Dphi, area, _) = gradbasis(mesh.p, mesh.t)



    NT = size(mesh.t, 1)
    Mt = zeros(NT, 3, 3)
    At = zeros(NT, 3, 3)
    for i = 1:3
        for j = 1:3        
            At[:, i, j] = (Dphi[:, 1, i].*Dphi[:, 1, j] + Dphi[:, 2, i].*Dphi[:, 2, j]).*area
            Mt[:, i, j] = area*((i==j)+1)/12
        end    
    end

    #%% Assemble the mass matrix in Omega
    M = spzeros(N, N)
    A = spzeros(N, N)
    for i = 1:3
        krow = mesh.t[:, i]
        for j = 1:3
            kcol = mesh.t[:, j]
            M = M + sparse(Int64.(krow), Int64.(kcol), Mt[:, i, j], N, N)
            A = A + sparse(Int64.(krow), Int64.(kcol), At[:, i, j], N, N)
        end 
    end
    return A, M, FreeNodes, dirichlet
end
#=
a=rect_grid2(1, 2, 1, 2, 10, 10)
Stiff_Mass(a)
=#

function gradbasis(node, elem)
    NT = size(elem, 1)
    elem = Int64.(elem)

    ve1 = node[elem[:, 3], :]-node[elem[:, 2], :]
    ve2 = node[elem[:, 1], :]-node[elem[:, 3], :]
    ve3 = node[elem[:, 2], :]-node[elem[:, 1], :]
    area = 0.5*(-ve3[:, 1].*ve2[:, 2] + ve3[:, 2].*ve2[:, 1])
    Dlambda = zeros(NT, 2, 3)
    Dlambda[1:NT, :, 3] = [-ve3[:, 2]./(2*area) ve3[:, 1]./(2*area)]
    Dlambda[1:NT, :, 1] = [-ve1[:, 2]./(2*area) ve1[:, 1]./(2*area)]
    Dlambda[1:NT, :, 2] = [-ve2[:, 2]./(2*area) ve2[:, 1]./(2*area)]

    idx = findall(x->x<0, area)
    #idx = (area<0)
    area[idx, :] = -area[idx, :]
    elemSign = ones(NT, 1)
    elemSign[idx] .= -1
    return Dlambda, area, elemSign
end

function f(p)        
    f = (sin.(2*pi*p[:, 1]).*sin.(2*pi*p[:, 2]))
    return f
end
function exactu(p, s)        
    u = (8*pi^2)^(-s)*sin.(2*pi*p[:, 1]).*sin.(2*pi*p[:, 2])
    return u
end
#s=0.2

function solve()

    for ind = 1:5
        mesh = rect_grid2(0, 1, 0, 1, 2^(ind+1), 2^(ind+1))

        (A, M, FreeNodes, dirichlet) = Stiff_Mass(mesh)

        N = size(A, 1)

        rhs = M*f(mesh.p)
        
        global u = zeros(N, 1)
        
        h = N^(-1/2)    

        k = 1/log(1/h)
        
        Nplus  = ceil(pi^2/(4*s*k^2))
        Nminus = ceil(pi^2/(4*(1-s)*k^2))
        
        for ell=-Nminus:1:Nplus
            y = k*ell;
            u[FreeNodes] = u[FreeNodes] + exp((1-s)*y) * ( (exp(y)*M[FreeNodes, FreeNodes]+A[FreeNodes, FreeNodes]) \ rhs[FreeNodes])
        end    
            
        u[FreeNodes] = (sin(s*pi)/pi) * k * u[FreeNodes]
        

        uD = zeros(N, 1)
        dirichlet = Int.(dirichlet)
        uD[dirichlet] = exactu(mesh.p[dirichlet, :], s)
        v = uD
        rhs = -A*uD
        v[FreeNodes] = A[FreeNodes, FreeNodes] \ rhs[FreeNodes]
        
        u = u+v
    
        #errL2[ind, 1] = sqrt((u-exactu(mesh.p, s))'*M*(u-exactu(mesh.p, s)))
        
        #hmax[ind, 1] = h
    end
    return u
end
#=
m = rect_grid2(0, 1, 0, 1, 2^(5+1), 2^(5+1))


using Plots
gr()
#pyplot()
#plotlyjs()
#pyplot()
plot3d(m.p[:, 1], m.p[:, 2], u)
=#